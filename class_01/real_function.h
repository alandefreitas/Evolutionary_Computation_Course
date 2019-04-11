//
// Created by Alan de Freitas on 05/04/2018.
//

#ifndef EVOLUTIONARY_COMPUTATION_REAL_FUNCTION_H
#define EVOLUTIONARY_COMPUTATION_REAL_FUNCTION_H

#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <algorithm>

class real_function {
public:
    // Evolutionary Algorithm
    // initialize with a rastrigin functions
    explicit real_function(size_t n);
    // initialize with any function and the same bound for all variable
    real_function(size_t n, std::function<double(const std::vector<double>&)> fx, double lower_bound, double upper_bound);
    // initialize with any function and bounds for each variable
    real_function(size_t n, std::function<double(const std::vector<double>&)> fx, std::vector<double> lower_bounds, std::vector<double> upper_bounds);
    void disp();
    size_t size();
    bool is_minimization() { return true; };
    // Auxiliary
    double lower_bound(size_t i);
    double upper_bound(size_t i);
    double fx(const std::vector<double>& x);
private:
    size_t _size;
    std::function<double(const std::vector<double>&)> fx_;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;
};


#endif //EVOLUTIONARY_COMPUTATION_REAL_FUNCTION_H
