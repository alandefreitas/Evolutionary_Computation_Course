//
// Created by Alan de Freitas on 12/04/2018.
//

#include "EA.h"

template <typename problem, typename solution>
std::default_random_engine EA<problem,solution>::_generator = std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count());

template <typename problem, typename solution>
EA<problem,solution>::EA(problem &p) :
    _problem(p),
    _mutation_strength(1.0/p.size()),
    _best_fx(p.is_minimization() ? std::numeric_limits<double>::max() : std::numeric_limits<double>::min())
{
    initialize_population();
}

template <typename problem, typename solution>
void EA<problem,solution>::initialize_population(){
    this->_population.reserve(_population_size);
    for (int i = 0; i < this->_population_size; ++i) {
        this->_population.emplace_back(std::make_shared<solution>(this->_problem));
    }
}

template <typename problem, typename solution>
void EA<problem,solution>::run() {
    for (int i = 0; i < this->_max_generations; ++i) {
        evolutionary_cycle();
    }
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
    evaluation_step();
    population_type children = reproduction_step();
    population_update_step(children);
}

template <typename problem, typename solution>
double EA<problem,solution>::best_fx() {
    return this->_best_fx;
}

template <typename problem, typename solution>
void EA<problem,solution>::evaluate(population_type& population){
    for (solution_ptr& item : population) {
        item->fx = item->evaluate(this->_problem);
        if (this->_problem.is_minimization()){
            if (item->fx < this->_best_fx){
                this->_best_fx = item->fx;
            }
        } else {
            if (item->fx > this->_best_fx){
                this->_best_fx = item->fx;
            }
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
typename EA<problem,solution>::population_type EA<problem,solution>::selection(population_type& population, size_t n_of_candidates, selection_strategy s){
    switch (s){
        case selection_strategy::uniform:
            return uniform_selection(population, n_of_candidates);
        case selection_strategy::truncate:
            return truncate_selection(population, n_of_candidates);
        case selection_strategy::tournament:
            return tournament_selection(population, n_of_candidates);
        case selection_strategy::roulette:
            return roulette_selection(population, n_of_candidates);
        case selection_strategy::sus:
            return sus_selection(population, n_of_candidates);
        case selection_strategy::overselection:
            return overselection_selection(population, n_of_candidates);
        case selection_strategy::roundrobin_tournament:
            return roundrobin_selection(population, n_of_candidates);
    }
};

template <typename problem, typename solution>
typename EA<problem,solution>::population_type EA<problem,solution>::reproduction(population_type& parents){
    std::uniform_real_distribution<double> r(0.0,1.0);
    population_type children;
    children.reserve(parents.size()/this->_parents_per_children);
    for (int j = 0; j < parents.size(); j += this->_parents_per_children) {
        if (r(EA::_generator) < this->_crossover_probability) {
            children.emplace_back(
                    std::make_shared<solution>(
                        parents[j]->crossover(
                            this->_problem,
                            *parents[j+1]
                        )
                    )
            );
        } else {
            children.emplace_back(std::make_shared<solution>(*parents[j]));
            children.back()->mutation(this->_problem,this->_mutation_strength);
        }
    }
    return children;
}

template <typename problem, typename solution>
void EA<problem,solution>::display_status() {
    std::cout << "Generation #" << ++_current_generation;
    if (this->_current_generation > 1){
        std::cout << " Best_fx: " << this->best_fx() << std::endl;
        double min_fx = std::numeric_limits<double>::max();
        double avg_fx = 0.0;
        double max_fx = -std::numeric_limits<double>::max();
        for (solution_ptr& ind : this->_population) {
            min_fx = std::min(min_fx,ind->fx);
            max_fx = std::max(max_fx,ind->fx);
            avg_fx += ind->fx;
        }
        avg_fx /= this->_population.size();
        std::cout << "| Min: " << min_fx
                  << " | Avg: " << avg_fx
                  << " | Max: " << max_fx
                  << " | "
                  << ((min_fx == max_fx) ? ("(Population Converged)"):(""))
                  << std::endl;
    } else {
        std::cout << std::endl;
    }
}


template <typename problem, typename solution>
void EA<problem,solution>::scaling(population_type& population, scaling_strategy s) {
    switch (s){
        case scaling_strategy::window: {
            window_scaling(population);
            break;
        }
        case scaling_strategy::sigma: {
            sigma_scaling(population);
            break;
        }
        case scaling_strategy::linear_rank: {
            linear_rank_scaling(population);
            break;
        }
        case scaling_strategy::exponential_rank: {
            exponential_rank_scaling(population);
            break;
        }
    }
}

template <typename problem, typename solution>
void EA<problem,solution>::attribute_fitness_from_rank(EA::population_type &population) {
    std::sort(population.begin(),population.end(),
              [](solution_ptr& a, std::shared_ptr<solution >&b){
                  return a->fx < b->fx;
              });
    for (int i = 0; i < population.size(); ++i) {
        population[i]->fitness = this->_problem.is_minimization() ? population.size() - i : i + 1;
    }
}

template <typename problem, typename solution>
typename EA<problem,solution>::population_type EA<problem,solution>::uniform_selection(EA::population_type &population, size_t n_of_candidates) {
    population_type parents(n_of_candidates);
    std::uniform_int_distribution<size_t> pos_d(0,population.size()-1);
    for (solution_ptr& parent : parents) {
        parent = population[pos_d(EA::_generator)];
    }
    return parents;

}

template <typename problem, typename solution>
typename EA<problem,solution>::population_type EA<problem,solution>::tournament_selection(EA::population_type &population, size_t n_of_candidates) {
    population_type parents;
    std::uniform_int_distribution<size_t> pos_d(0,population.size()-1);
    for (int i = 0; i < n_of_candidates; ++i) {
        size_t position = pos_d(EA::_generator);
        for (int j = 1; j < this->_tournament_size; ++j) {
            size_t position2 = pos_d(EA::_generator);
            if (population[position2]->fitness > population[position]->fitness) {
                position = position2;
            }
        }
        parents.push_back(population[position]);
    }
    return parents;
}

template <typename problem, typename solution>
typename EA<problem,solution>::population_type EA<problem,solution>::roulette_selection(EA::population_type &population, size_t n_of_candidates) {
    population_type parents;
    parents.reserve(n_of_candidates);
    std::discrete_distribution<size_t> pos_d (population.size(),0,population.size()-1,[&population](size_t pos){
        return population[pos]->fitness;
    });
    for (int i = 0; i < n_of_candidates; ++i) {
        parents.push_back(population[pos_d(EA::_generator)]);
    }
    return parents;
}

template <typename problem, typename solution>
typename EA<problem,solution>::population_type EA<problem,solution>::truncate_selection(EA::population_type &population, size_t n_of_candidates) {
    population_type parents;
    parents.reserve(n_of_candidates);
    std::partial_sort(population.begin(),
                      population.begin() + std::min(n_of_candidates,population.size()),
                      population.end(),
                      [](solution_ptr& a, solution_ptr& b){
                          return a->fitness > b->fitness;
                      }
    );
    std::copy(population.begin(),population.begin() + std::min(n_of_candidates,population.size()),std::back_inserter(parents));
    int i = 0;
    while (parents.size() < n_of_candidates){
        parents.push_back(parents[i%population.size()]);
        i++;
    }
    std::shuffle(parents.begin(),parents.end(),EA::_generator);
    return parents;
}

template <typename problem, typename solution>
typename EA<problem,solution>::population_type EA<problem,solution>::sus_selection(EA::population_type &population, size_t n_of_candidates) {
    population_type parents;
    parents.reserve(n_of_candidates);
    double total_fit = 0.0;
    for (solution_ptr& ind : population) {
        total_fit += ind->fitness;
    }
    double gap = total_fit/n_of_candidates;
    std::uniform_real_distribution<double> dist_r(0.0,gap);
    double r = dist_r(EA::_generator);
    size_t current_ind = 0;
    double sum = population[current_ind]->fitness;
    for (int i = 0; i < n_of_candidates; ++i) {
        while (r > sum){
            ++current_ind;
            sum += population[current_ind]->fitness;
        }
        parents.push_back(population[current_ind]);
        r += gap;
    }
    std::shuffle(parents.begin(),parents.end(),EA::_generator);
    return parents;
}

template <typename problem, typename solution>
typename EA<problem,solution>::population_type EA<problem,solution>::overselection_selection(EA::population_type &population, size_t n_of_candidates) {
    population_type parents;
    parents.reserve(n_of_candidates);
    std::partial_sort(population.begin(),
                      population.begin() + population.size()*this->_overselection_proportion,
                      population.end(),
                      [](solution_ptr& a, solution_ptr& b){
                          return a->fitness > b->fitness;
                      }
    );
    std::uniform_int_distribution<size_t> pos_best(0,population.size()*this->_overselection_proportion-1);
    for (int i = 0; i < n_of_candidates * 0.8; ++i) {
        parents.push_back(population[pos_best(EA::_generator)]);
    }
    std::uniform_int_distribution<size_t> pos_worst(population.size()*this->_overselection_proportion,population.size()-1);
    for (int i = 0; i < n_of_candidates - (n_of_candidates * 0.8); ++i) {
        parents.push_back(population[pos_worst(EA::_generator)]);
    }
    return parents;
}

template <typename problem, typename solution>
typename EA<problem,solution>::population_type EA<problem,solution>::roundrobin_selection(EA::population_type &population, size_t n_of_candidates) {
    struct solution_with_score {
        solution_ptr s;
        size_t score = 0;
    };
    std::vector<solution_with_score> scores(population.size());
    for (int i = 0; i < population.size(); ++i) {
        scores[i].s = population[i];
    }
    std::uniform_int_distribution<size_t> pos_d(0,scores.size()-1);
    for (int i = 0; i < population.size(); ++i) {
        solution_with_score& player1 = scores[i];
        for (int j = 1; j < this->_roundrobin_tournament_size; ++j) {
            solution_with_score& player2 = scores[pos_d(EA::_generator)];
            if (player1.s->fitness > player2.s->fitness) {
                player1.score++;
            } else {
                player2.score++;
            }
        }
    }
    std::sort(scores.begin(),
              scores.end(),
              [](solution_with_score& a, solution_with_score& b){
                  return a.score > b.score;
              }
    );
    population_type parents;
    parents.reserve(n_of_candidates);
    for (int i = 0; i < std::min(n_of_candidates,population.size()); ++i) {
        parents.push_back(scores[i].s);
    }
    int i = 0;
    while (parents.size() < n_of_candidates){
        parents.push_back(parents[i%population.size()]);
        i++;
    }
    return parents;
}

template <typename problem, typename solution>
void EA<problem,solution>::window_scaling(EA::population_type &population) {
    for (solution_ptr& ind : population) {
        ind->fitness = this->_problem.is_minimization() ? -ind->fx : ind->fx;
    }
    const auto iter_to_min_fitness = std::min_element(population.begin(),
                                                      population.end(),
                                                      [](solution_ptr& a,
                                                         std::shared_ptr<solution >&b){
                                                          return a->fitness < b->fitness;
                                                      });
    double min_fitness = (*iter_to_min_fitness)->fitness;
    for (solution_ptr& ind : population) {
        ind->fitness -= min_fitness + 1;
    }
}

template <typename problem, typename solution>
void EA<problem,solution>::sigma_scaling(EA::population_type &population) {
    for (solution_ptr& ind : population) {
        ind->fitness = this->_problem.is_minimization() ? -ind->fx : ind->fx;
    }
    double avg_fitness = 0.0;
    for (solution_ptr& ind : population) {
        avg_fitness += ind->fitness;
    }
    double std_fitness = 0.0;
    for (solution_ptr& ind : population) {
        std_fitness += pow(ind->fitness - avg_fitness,2.0);
    }
    std_fitness /= population.size() - 1;
    std_fitness = sqrt(std_fitness);
    for (solution_ptr& ind : population) {
        ind->fitness = std::max(
                this->_sigma_bias + (ind->fitness - avg_fitness)/(this->_sigma_constant * std_fitness),
                0.0);
    }
}

template <typename problem, typename solution>
void EA<problem,solution>::linear_rank_scaling(EA::population_type &population) {
    attribute_fitness_from_rank(population);
    const double bias = ((2-_linear_rank_selective_pressure)/population.size());
    const double denominator = population.size() * (population.size() - 1);
    for (solution_ptr& ind : population) {
        ind->fitness = bias + (2 * ind->fitness * (_linear_rank_selective_pressure - 1))/denominator;
    }
}

template <typename problem, typename solution>
void EA<problem,solution>::exponential_rank_scaling(EA::population_type &population) {
    attribute_fitness_from_rank(population);
    for (solution_ptr& ind : population) {
        ind->fitness = (1-exp(-ind->fitness))/population.size();
    }
}

template <typename problem, typename solution>
void EA<problem,solution>::evaluation_step() {
    evaluate(this->_population);
    scaling(this->_population,scaling_strategy::exponential_rank);
}

template <typename problem, typename solution>
typename EA<problem,solution>::population_type EA<problem,solution>::reproduction_step() {
    population_type parents = selection(this->_population, n_of_selection_candidates(),selection_strategy::roundrobin_tournament);
    population_type children = reproduction(parents);
    evaluate(children);
    return children;
}

template <typename problem, typename solution>
void EA<problem,solution>::population_update_step(population_type& children) {
    // update population
    population_type next_population_candidates = children;
    if (this->_competition_between_parents_and_children ||
        this->_children_proportion < 1.0){
        next_population_candidates.insert(
                next_population_candidates.end(),
                this->_population.begin(),
                this->_population.end()
        );
    }
    size_t size_of_elite_set = size_t(ceil(this->_elitism_proportion * this->_population.size()));
    scaling(next_population_candidates,scaling_strategy::window);
    next_population_candidates = selection(next_population_candidates, this->_population_size - size_of_elite_set,selection_strategy::truncate);
    // add the elite set to the next population
    std::partial_sort(this->_population.begin(),
                      this->_population.begin() + size_of_elite_set,
                      this->_population.end(),
                      [](solution_ptr& a, solution_ptr& b){
                          return a->fitness > b->fitness;
                      }
    );
    next_population_candidates.insert(
            next_population_candidates.end(),
            this->_population.begin(),
            this->_population.begin() + size_of_elite_set
    );
    this->_population = next_population_candidates;
}
