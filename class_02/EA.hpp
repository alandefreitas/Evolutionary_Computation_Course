//
// Created by Alan de Freitas on 12/04/2018.
//

#include "EA.h"

template <typename problem, typename solution>
std::default_random_engine EA<problem,solution>::_generator = std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count());

template <typename problem, typename solution>
EA<problem,solution>::EA(problem &p) :
    _problem(p),
    // Parameters
    _max_generations(1000),
    _population_size(200),
    _parents_per_children(2),
    _children_proportion(7),
    _crossover_probability(0.9),
    _mutation_strength(1.0/this->_problem.size()),
    // Data
    _current_generation(0)
{
    // Initialize population
    this->_population.reserve(_population_size);
    for (int i = 0; i < this->_population_size; ++i) {
        this->_population.emplace_back(this->_problem);
    }
}

template <typename problem, typename solution>
void EA<problem,solution>::run() {
    for (int i = 0; i < this->_max_generations; ++i) {
        evolutionary_cycle();
    }
    // To be replaced by:
    //     while (!stopping_criteria()){
    //         evolutionary_cycle();
    //     }
}

template <typename problem, typename solution>
void EA<problem,solution>::run(size_t iterations) {
    for (int i = 0; i < iterations; ++i) {
        evolutionary_cycle();
    }
}

template <typename problem, typename solution>
void EA<problem,solution>::evolutionary_cycle() {
    display_status();
    evaluate(this->_population);
    std::vector<size_t> parent_position = selection(this->_population, n_of_selection_candidates(),selection_strategy::uniform);
    std::vector<solution> children = reproduction(this->_population, parent_position);
    evaluate(children);
    std::vector<size_t> children_position = selection(this->_population, this->_population_size,selection_strategy::truncate);
    this->_population = update_population(children,children_position);
}


template <typename problem, typename solution>
double EA<problem,solution>::best_fx() {
    if (this->_problem.is_minimization()) {
        return -this->_best_fx;
    } else {
        return this->_best_fx;
    }
}

template <typename problem, typename solution>
void EA<problem,solution>::evaluate(std::vector<solution>& population){
    for (solution& item : population) {
        item.fx = item.evaluate(this->_problem);
        if (this->_problem.is_minimization()){
            item.fx = -item.fx;
        }
        if (item.fx > this->_best_fx){
            this->_best_fx = item.fx;
        }
    }
};

template <typename problem, typename solution>
size_t EA<problem,solution>::n_of_selection_candidates(){
    return this->_population_size*
           this->_parents_per_children*
           this->_children_proportion;
};

template <typename problem, typename solution>
std::vector<size_t> EA<problem,solution>::selection(std::vector<solution>& population,
                                                    size_t n_of_candidates,
                                                    selection_strategy s){
    switch (s){
        case selection_strategy::uniform: {
            std::vector<size_t> parent_position(n_of_candidates);
            std::uniform_int_distribution<size_t> pos_d(0,population.size()-1);
            for (size_t& position : parent_position) {
                position = pos_d(EA::_generator);
            }
            return parent_position;
        }
        case selection_strategy::truncate: {
            std::vector<size_t> parent_position(n_of_candidates);
            std::partial_sort(population.begin(),
                              population.begin() + parent_position.size(),
                              population.end(),
                              [](solution& a, solution& b){
                                  return a.fx > b.fx;
                              }
            );
            std::iota(parent_position.begin(),parent_position.end(),0);
            return parent_position;
        }
    }
};

template <typename problem, typename solution>
std::vector<solution> EA<problem,solution>::reproduction(std::vector<solution>& population, std::vector<size_t>& parent_position){
    std::uniform_real_distribution<double> r(0.0,1.0);
    std::vector<solution> children;
    for (int j = 0; j < parent_position.size(); j += 2) {
        if (r(EA::_generator) < this->_crossover_probability) {
            // Crossover
            children.push_back(
                    population[parent_position[j]].crossover(
                            this->_problem,
                            population[parent_position[j+1]]
                    )
            );
        } else {
            // Mutation
            children.push_back(population[parent_position[j]]);
            children.back().mutation(this->_problem,this->_mutation_strength);
        }
    }
    return children;
}

template <typename problem, typename solution>
std::vector<solution> EA<problem,solution>::update_population(std::vector<solution>& population, std::vector<size_t>& positions) {
    std::vector<solution> r;
    r.reserve(population.size());
    for (size_t position : positions) {
        r.push_back(population[position]);
    }
    return r;
}

template <typename problem, typename solution>
void EA<problem,solution>::display_status() {
    std::cout << "Generation #" << ++_current_generation;
    std::cout << " - Best_fx: " << this->best_fx() << std::endl;
}
