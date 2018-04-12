//
// Created by Alan de Freitas on 05/04/2018.
//

#ifndef EVOLUTIONARY_COMPUTATION_KNAPSACKP_H
#define EVOLUTIONARY_COMPUTATION_KNAPSACKP_H

#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <algorithm>

class knapsack_p {
public:
    // Evolutionary Algorithm
    knapsack_p(size_t n);
    void disp();
    size_t size();
    // Problem-specific
    int value(size_t i);
    int weight(size_t i);
    int capacity();
    int max_value();
private:
    int _capacity;
    int _max_value;
    std::vector<int> _value;
    std::vector<int> _weight;
    static std::default_random_engine _generator;
};


#endif //EVOLUTIONARY_COMPUTATION_KNAPSACKP_H
