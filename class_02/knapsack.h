//
// Created by Alan de Freitas on 05/04/2018.
//

#ifndef EVOLUTIONARY_COMPUTATION_KNAPSACK_H
#define EVOLUTIONARY_COMPUTATION_KNAPSACK_H

#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <algorithm>

#include "knapsackP.h"

class knapsack {
public:
    knapsack(knapsack_p &p);
    void disp(knapsack_p &p);
    double evaluate(knapsack_p &p);
    void mutation(knapsack_p &p, double mutation_strength);
    knapsack crossover(knapsack_p &p, knapsack& rhs);
    double fx;
private:
    std::vector<int> _knapsack;
    static std::default_random_engine _generator;
};

#endif //EVOLUTIONARY_COMPUTATION_KNAPSACK_H
