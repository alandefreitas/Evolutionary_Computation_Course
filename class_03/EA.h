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
    enum selection_strategy {uniform, truncate, tournament, roulete, sus};
    enum scaling_strategy {sigma, rank};
public:
    EA(problem &p);
    void run();
    void run(size_t iterations);
    double best_fx();
private:
    void evolutionary_cycle();
    void evaluate(std::vector<solution>& population);
    void scaling(std::vector<solution>& population, scaling_strategy s);
    size_t n_of_selection_candidates();
    std::vector<size_t> selection(std::vector<solution>& population, size_t n_of_candidates, selection_strategy s);
    std::vector<solution> reproduction(std::vector<solution>& population, std::vector<size_t>& parent_position);
    std::vector<solution> update_population(std::vector<solution>& population, std::vector<size_t>& positions);
    void display_status();

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
    size_t _current_generation;

    // Auxiliary generator
    static std::default_random_engine _generator;

    // Get statistics
    double _best_fx = -std::numeric_limits<double>::max();

};

#endif //EVOLUTIONARY_COMPUTATION_EA_H
