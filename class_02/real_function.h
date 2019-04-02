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
    real_function(size_t n);
    void disp();
    size_t size();
    bool is_minimization() { return true; };
    // Auxiliary
    double lower_bound() { return -5.12; }
    double upper_bound() { return +5.12; }
private:
    size_t _size;
};


#endif //EVOLUTIONARY_COMPUTATION_REAL_FUNCTION_H
