//
// Created by Alan de Freitas on 05/04/2018.
//

#ifndef EVOLUTIONARY_COMPUTATION_REAL_VECTOR_H
#define EVOLUTIONARY_COMPUTATION_REAL_VECTOR_H

#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <algorithm>

#include "real_function.h"

class real_vector {
public:
    real_vector(real_function &p);
    void disp(real_function &p);
    double evaluate(real_function &p);
    void mutation(real_function &p, double mutation_strength);
    real_vector crossover(real_function &p, real_vector& rhs);
    double fx;
private:
    std::vector<double> _real_vector;
    static std::default_random_engine _generator;
};

#endif //EVOLUTIONARY_COMPUTATION_REAL_VECTOR_H
