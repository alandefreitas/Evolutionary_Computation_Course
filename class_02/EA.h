//
// Created by Alan de Freitas on 12/04/2018.
//

#ifndef EVOLUTIONARY_COMPUTATION_EA_H
#define EVOLUTIONARY_COMPUTATION_EA_H

#include "knapsackP.h"
#include "knapsack.h"

template <typename problem, typename solution>
class EA {
public:
    EA(problem &p);
    void run();
    double max_fx();
private:
    // Search parameters
    size_t _max_generations;
    size_t _population_size;
    size_t _parents_per_children;
    size_t _children_proportion;
    double _crossover_probability;
    double _mutation_strength;

    // Data
    std::vector<solution> _population;
    problem _problem;

    // Auxiliary generator
    static std::default_random_engine _generator;

    // Get statistics
    double _max_fx = std::numeric_limits<double>::min();
};

#endif //EVOLUTIONARY_COMPUTATION_EA_H
