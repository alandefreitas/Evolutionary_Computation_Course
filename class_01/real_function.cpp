//
// Created by Alan de Freitas on 05/04/2018.
//

#include "real_function.h"


real_function::real_function(size_t n)
        : _size(n),
        lower_bounds_(n,-5.12),
        upper_bounds_(n,+5.12)
{
    // Since the user gave no parameters, we initialize it with a rastrigin function
    fx_ = [n](const std::vector<double> &x) {
        const double A = 10;
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
            sum += pow(x[i],2) - A * cos(2. * 3.14 * x[i]);
        }
        return A * n + sum;
    };
}

real_function::real_function(size_t n, std::function<double(const std::vector<double>&)> fx, double lower_bound, double upper_bound)
        : _size(n),
        fx_(fx),
        lower_bounds_(n,lower_bound),
        upper_bounds_(n,upper_bound)
{ }

real_function::real_function(size_t n, std::function<double(const std::vector<double>&)> fx, std::vector<double> lower_bounds, std::vector<double> upper_bounds)
        : _size(n),
        fx_(fx),
        lower_bounds_(lower_bounds),
        upper_bounds_(upper_bounds)
{}

void real_function::disp() {
    std::cout << "Real function of size " << _size << std::endl;
}

size_t real_function::size() {
    return this->_size;
}

double real_function::lower_bound(size_t i) {
    return lower_bounds_[i];
}

double real_function::upper_bound(size_t i) {
    return upper_bounds_[i];
}

double real_function::fx(const std::vector<double>& x){
    return this->fx_(x);
}