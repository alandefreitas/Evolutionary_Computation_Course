//
// Created by Alan de Freitas on 12/04/2018.
//

#ifndef EVOLUTIONARY_COMPUTATION_EA_H
#define EVOLUTIONARY_COMPUTATION_EA_H

#include "knapsackP.h"
#include "knapsack.h"

class EA {
public:
    EA(knapsack_p &p);
    void run();
    double max_fx();
private:
    // Search parameters
    size_t max_generations;
    size_t population_size;
    size_t parents_per_children;
    size_t children_proportion;
    double crossover_probability;
    double mutation_strength;

    // Data
    std::vector<knapsack> population;
    knapsack_p problem;

    // Auxiliary generator
    static std::default_random_engine _generator;

    // Get statistics
    double max = std::numeric_limits<double>::min();
};


#endif //EVOLUTIONARY_COMPUTATION_EA_H
