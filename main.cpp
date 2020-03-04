#include <iostream>
#include <vector>
#include "integrals.h"
#include <limits>
#include <cmath>

typedef std::numeric_limits< double > dbl;

struct run_params_t {
    double e_r;
    double e_a;
    int n_threads;
    double x_start; double x_end;
    double y_start; double y_end;
};

void benchmark(run_params_t* arguments);

int main() {
    std::cout.precision(dbl::digits10);
    run_params_t pars{pow(10, -4), pow(10, 1),
                      7, -101, 101, -101 ,101};
    benchmark(&pars);
    return 0;
}


void benchmark(run_params_t* arguments) {
    auto ars = *arguments;
    double local_e_r = dbl::infinity();
    double local_e_a = dbl::infinity();
    int steps = 1171;

    integral_data_t params
    {
    ars.x_start, ars.x_end,
    ars.y_start, ars.y_end,
    steps
    };

    if (arguments-> n_threads == 1) {
        no_paralel(&params);
    }else {
        paralelize(&params, arguments->n_threads);
    }
    auto prev_res = params.res;

    while (local_e_a > arguments->e_a && local_e_r > arguments->e_r) {

        params.n_steps += 1000;
        if (arguments-> n_threads == 1) {
            no_paralel(&params);
        }else {
            paralelize(&params, arguments->n_threads);
        }

        local_e_a = fabs(params.res - prev_res);
        local_e_r = local_e_a / prev_res;
        prev_res = params.res;
    }
    std::cout << "Result = " << params.res << std::endl;
    std::cout << "E_rel = " << local_e_r << "\tE_abs = " << local_e_a << std::endl;
    std::cout << "Time = " << params.time << "ms" << std::endl;
}