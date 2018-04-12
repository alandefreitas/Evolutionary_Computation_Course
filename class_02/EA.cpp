//
// Created by Alan de Freitas on 12/04/2018.
//

#include "EA.h"

std::default_random_engine EA::_generator = std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count());

EA::EA(knapsack_p &p) :
    problem(p) {
    // Parameters
    max_generations = 1000;
    population_size = 200;
    parents_per_children = 2;
    children_proportion = 7;
    crossover_probability = 0.9;
    mutation_strength = 1/problem.size();

    // Initialize population
    population.reserve(population_size);
    for (int i = 0; i < population_size; ++i) {
        population.emplace_back(problem);
    }
}

void EA::run() {
    // Evolutionary Cycle
    for (int i = 0; i < max_generations; ++i) {
        std::cout << "Generation #" << i + 1 << " fx: " << this->max_fx() << std::endl;
        // Evaluate
        for (knapsack& item : population) {
            item.fx = item.evaluate(problem);
        }
        // Selection
        std::vector<size_t> parent_position(population_size*parents_per_children*children_proportion);
        std::uniform_int_distribution<size_t> pos_d(0,population_size-1);
        for (size_t& position : parent_position) {
            position = pos_d(EA::_generator);
        }
        // Reproduction
        std::uniform_real_distribution<double> r(0,1);
        std::vector<knapsack> children;
        for (int j = 0; j < parent_position.size(); j += 2) {
            if (pos_d(EA::_generator) < crossover_probability) {
                // Crossover
                children.push_back(
                        population[parent_position[j]].crossover(
                                problem,population[parent_position[j+1]]
                        )
                );
            } else {
                // Mutation
                children.push_back(population[parent_position[j]]);
                children.back().mutation(problem,mutation_strength);
            }
        }
        // Survival
        std::partial_sort(children.begin(),
                          children.begin()+population_size,
                          children.end(),
                          [](knapsack& a, knapsack& b){
                              return a.fx > b.fx;
                          }
        );
        // Save best fx
        max = std::max(max,children.front().fx);
    }
}

double EA::max_fx() {
    return this->max;
}
