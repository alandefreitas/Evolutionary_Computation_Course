//
// Created by Alan de Freitas on 12/04/2018.
//

#ifndef EVOLUTIONARY_COMPUTATION_NAIVE_EA_H
#define EVOLUTIONARY_COMPUTATION_NAIVE_EA_H

#include "knapsackP.h"
#include "knapsack.h"

template<typename problem, typename solution>
class naive_EA {
public:
    enum selection_strategy {
        uniform, truncate
    };
public:
    // Auxiliary class for individuals
    class individual {
    public:
        individual(solution rhs) : s(rhs){}
        individual(problem rhs) : s(rhs){}
        double evaluate(problem& rhs) { return s.evaluate(rhs); }
        solution s;
        double fx;
    };
public:
    explicit naive_EA(problem &p);

    void run();

    void run(size_t iterations);

    double best_fx();

private:
    void evolutionary_cycle();

    void evaluate(std::vector<individual> &population);

    size_t n_of_selection_candidates();

    std::vector<size_t> selection(std::vector<individual> &population, size_t n_of_candidates, selection_strategy s);

    std::vector<individual> reproduction(std::vector<individual> &population, std::vector<size_t> &parent_position);

    std::vector<individual> update_population(std::vector<individual> &population, std::vector<size_t> &positions);

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
    std::vector<individual> _population;
    problem _problem;
    size_t _current_generation;

    // Auxiliary generator
    static std::default_random_engine _generator;

    // Get statistics
    double _best_fx = -std::numeric_limits<double>::max();
};

#endif //EVOLUTIONARY_COMPUTATION_NAIVE_EA_H
