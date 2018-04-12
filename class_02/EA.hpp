//
// Created by Alan de Freitas on 12/04/2018.
//

#include "EA.h"

template <typename problem, typename solution>
std::default_random_engine EA<problem,solution>::_generator = std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count());

template <typename problem, typename solution>
EA<problem,solution>::EA(problem &p) :
    _problem(p) {
    // Parameters
    this->_max_generations = 1000;
    this->_population_size = 200;
    this->_parents_per_children = 2;
    this->_children_proportion = 7;
    this->_crossover_probability = 0.9;
    this->_mutation_strength = 1/this->_problem.size();

    // Initialize population
    this->_population.reserve(_population_size);
    for (int i = 0; i < this->_population_size; ++i) {
        this->_population.emplace_back(this->_problem);
    }
}

template <typename problem, typename solution>
void EA<problem,solution>::run() {
    // Evolutionary Cycle
    for (int i = 0; i < this->_max_generations; ++i) {
        std::cout << "Generation #" << i + 1 << " fx: " << this->max_fx() << std::endl;
        // Evaluate
        for (solution& item : this->_population) {
            item.fx = item.evaluate(this->_problem);
        }
        // Selection
        const size_t n_of_candidates = this->_population_size*
                this->_parents_per_children*
                this->_children_proportion;
        std::vector<size_t> parent_position(n_of_candidates);
        std::uniform_int_distribution<size_t> pos_d(0,this->_population_size-1);
        for (size_t& position : parent_position) {
            position = pos_d(EA::_generator);
        }
        // Reproduction
        std::uniform_real_distribution<double> r(0.0,1.0);
        std::vector<solution> children;
        for (int j = 0; j < parent_position.size(); j += 2) {
            if (r(EA::_generator) < this->_crossover_probability) {
                // Crossover
                children.push_back(
                        this->_population[parent_position[j]].crossover(
                                this->_problem,
                                this->_population[parent_position[j+1]]
                        )
                );
            } else {
                // Mutation
                children.push_back(this->_population[parent_position[j]]);
                children.back().mutation(this->_problem,this->_mutation_strength);
            }
        }
        // Evaluate
        for (solution& item : children) {
            item.fx = item.evaluate(this->_problem);
        }
        // Survival
        std::partial_sort(children.begin(),
                          children.begin()+this->_population_size,
                          children.end(),
                          [](solution& a, solution& b){
                              return a.fx > b.fx;
                          }
        );
        // Save best fx
        this->_max_fx = std::max(this->_max_fx,children.front().fx);
    }
}

template <typename problem, typename solution>
double EA<problem,solution>::max_fx() {
    return this->_max_fx;
}