//
// Created by vincent on 02.03.20.
//

#ifndef LAB3_INTEGRALS_H
#define LAB3_INTEGRALS_H

#include <vector>

struct integral_data_t {
    double x_start;
    double x_end;
    double y_start;
    double y_end;
    int n_steps;
    double res;
    long long time;
};

//std::vector<std::vector<double>> divide_range(double start, double end, double n_steps, int n_ranges) ;

std::vector<std::vector<double>> divide_range(double start, double end, int n_ranges);
//void calc_integral(double x_begin, double x_end, double y_begin, double y_end, int steps , double * res_ptr, long long * time_ptr);
void calc_integral(integral_data_t* arguments);
//double paralelize(std::vector<int>* run_params, std::vector<double>* coords );

void paralelize(integral_data_t* arguments, int n_threads);
double no_paralel(integral_data_t* arguments);
double f(double x1, double x2);


#endif //LAB3_INTEGRALS_H
