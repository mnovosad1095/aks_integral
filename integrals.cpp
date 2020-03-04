//
// Created by vincent on 02.03.20.
//

#include "integrals.h"
#include <vector>
#include <thread>
#include <iostream>
#include <cmath>


inline std::chrono::high_resolution_clock::time_point get_current_time_fenced()
{
    std::atomic_thread_fence(std::memory_order_seq_cst);
    auto res_time = std::chrono::high_resolution_clock::now();
    std::atomic_thread_fence(std::memory_order_seq_cst);
    return res_time;
}

template<class D>
inline long long to_mili(const D& d)
{
    return std::chrono::duration_cast<std::chrono::milliseconds>(d).count();
}

std::vector<std::vector<double>> divide_range(double start, double end, int n_ranges) {
    std::vector<std::vector<double>> ranges;
    std::vector<double> pivots;
    auto delta = (end - start)/ n_ranges;

    pivots.push_back(start);
    for (int i = 1; i <= n_ranges;++i) {
        pivots.push_back(pivots[i-1] + delta);
    }

    for (int i = 1; i <= n_ranges; ++i) {
        ranges.push_back(std::vector<double>{pivots[i-1], pivots[i] - pow(10, -6)});
    }

    ranges[n_ranges-1][1] = end;

    return ranges;
}


void paralelize(integral_data_t* arguments, int n_threads){
    auto n_steps_by_thread = arguments->n_steps / n_threads;
    double local_res = 0;
    long long local_time = 0;
    auto max_time = local_time;

    std::vector<std::thread> th;
    th.reserve(n_threads);
    std::vector<integral_data_t> data(n_threads);
    std::vector<std::vector<double>> x_ranges;
    std::vector<std::vector<double>> y_ranges;

    if (arguments->x_end - arguments->x_start <= arguments->y_end - arguments->y_start) {
        x_ranges = divide_range(arguments->x_start, arguments->x_end, 2);
        y_ranges = divide_range(arguments->y_start, arguments->y_end, n_threads / 2);
    } else {
        x_ranges = divide_range(arguments->x_start, arguments->x_end, n_threads / 2);
        y_ranges = divide_range(arguments->y_start, arguments->y_end, 2);
    }

    int i = 0;
    for (auto &x_r: x_ranges) {
        for (auto &y_r: y_ranges) {

            data[i] = {x_r[0], x_r[1],
                      y_r[0], y_r[1],
                      n_steps_by_thread};
            th.emplace_back(calc_integral, &data[i]);
            ++i;
        }
    }

    for (auto &t: th) {
        t.join();
    }

    if (arguments->res != 0.0) {
        arguments->res = 0;
    }

    for (int j = 0; j < n_threads ; ++j) {
        arguments->res += data[j].res;
//        local_res += data[j].res;
        if (data[j].time > max_time) {
            max_time = data[j].time;
        }
    }

//    arguments->res = local_res;
    arguments->time = max_time;
}

double no_paralel(integral_data_t* arguments) {
    calc_integral(arguments);
    return 0;
}

void calc_integral(integral_data_t* arguments) {
    double local_res = 0;
    auto delta_x = (arguments->x_end - arguments->x_start) / arguments->n_steps;
    auto delta_y = (arguments->y_end - arguments->y_start) / arguments->n_steps;
    auto x = arguments->x_start;
    auto y = arguments->y_start;


    auto start = get_current_time_fenced();
    for (int i = 0; i < arguments->n_steps; ++i) {
        for(int j = 0; j < arguments->n_steps; ++j) {
            //calculate volume in current position
            local_res += delta_x * delta_y * f(x,y);
            // increase y
            y += delta_y;
        }
        //increase x
        x += delta_x;
        y = arguments->y_start;
    }
    auto final = get_current_time_fenced();
    arguments->time = to_mili(final - start);
    arguments->res = local_res;
}

double f(double x1, double x2) {
    double a = 20;
    double b = 0.2;
    double c = 2*M_PI;
    return -a*exp(-b * sqrt((x1*x1+ x2*x2)/2)) -
           exp((cos(c*x1)+cos(c*x2))/2) + a + exp(1);
}

#pragma clang diagnostic pop