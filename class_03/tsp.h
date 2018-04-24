//
// Created by Alan de Freitas on 05/04/2018.
//

#ifndef EVOLUTIONARY_COMPUTATION_TSP_H
#define EVOLUTIONARY_COMPUTATION_TSP_H

#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <algorithm>

class tsp {
public:
    // Evolutionary Algorithm
    tsp(size_t n);
    void disp();
    size_t size();
    bool is_minimization() { return true; };
    // Problem-specific
    double distance(size_t i, size_t j);
private:
    std::vector<std::vector<double>> _distance;
    static std::default_random_engine _generator;
};


#endif //EVOLUTIONARY_COMPUTATION_TSP_H
