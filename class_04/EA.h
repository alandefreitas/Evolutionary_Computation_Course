//
// Created by Alan de Freitas on 12/04/2018.
//

#ifndef EVOLUTIONARY_COMPUTATION_EA_H
#define EVOLUTIONARY_COMPUTATION_EA_H

#include <memory>
#include <cmath>
#include <functional>

template<typename problem, typename solution>
class EA {
    public:
        using solution_ptr = std::shared_ptr<solution>;
        using population_type = std::vector<std::shared_ptr<solution>>;

        // how solutions can be selected for reproduction or survival
        enum class selection_strategy {
            uniform, // same probability to all solutions
            truncate, // only the solutions with best fitness
            tournament, // winners of tournaments between groups of solutions
            roulette, // selection probability proportional to fitness
            sus, // stochastic probability proportional to fitness
            overselection, // 80% of operations on the x% best solutions
            roundrobin_tournament // each individual evaluated against other q individuals
        };

        // how fitness is determined from the objective functions
        enum class scaling_strategy {
            window, // fitness <- fx - worst_fx
            sigma, // fitness <- max(bias + (fx - avg_fx)/(c * std_fx)), 0.0)
            linear_rank, // fitness <- ((2-pressure_constant)/n_ind) + ((2 * rank(fx) * (pressure_constant - 1))/n_ind * (n_ind - 1))
            exponential_rank // fitness <- (1-e^(-i))/n_ind
        };

        // the structure of the islands
        enum class island_structure {
            ring, // each island is connected to the next
            random // each island is connected to another random island
        };

        // which individuals migrate from islands
        enum class island_migration_policy {
            best, // the best individual moves to another island
            fittest_half, // the best individual moves to another island
            random // a random individual moves to another island
        };

        // which individuals migrate from islands
        enum class island_replacement_policy {
            worst_swap,
            worst_overwrite,
            random_swap,
            random_overwrite,
        };

    public:
        EA(problem &p);

        void run();

        void run(size_t iterations);

        double best_fx();

        solution_ptr best_solution();

        typename population_type::iterator best_solutions_begin();

        typename population_type::iterator best_solutions_end();

    private /* methods */:
        // main cycle
        void evolutionary_cycle();

        // helpers
        size_t n_of_selection_candidates();

        void display_status();

        // evaluation algorithms
        void evaluation_step();

        void evaluate(population_type &population);

        // scaling algorithms
        void scaling(population_type &population, scaling_strategy s);

        void window_scaling(population_type &population);

        void sigma_scaling(population_type &population);

        void linear_rank_scaling(population_type &population);

        void exponential_rank_scaling(population_type &population);

        // fitness sharing
        void fitness_sharing(population_type &population);

        // selection algorithms
        population_type selection(population_type &population, size_t n_of_candidates, selection_strategy s);

        population_type uniform_selection(population_type &population, size_t n_of_candidates);

        population_type truncate_selection(population_type &population, size_t n_of_candidates);

        population_type tournament_selection(population_type &population, size_t n_of_candidates);

        population_type roulette_selection(population_type &population, size_t n_of_candidates);

        population_type sus_selection(population_type &population, size_t n_of_candidates);

        population_type overselection_selection(population_type &population, size_t n_of_candidates);

        population_type roundrobin_selection(population_type &population, size_t n_of_candidates);

        // reproduction algorithms
        population_type reproduction_step();

        population_type reproduction(population_type &parents);

        // update
        void population_update_step(population_type &children);

        void initialize_population();

        void attribute_fitness_from_rank(population_type &population);

        population_type get_island(population_type &population, int idx);

        void try_to_update_best(solution_ptr &candidate);

        void migration_step();

    private /* members */:
        // Search parameters
        // Population management
        size_t _population_size = 300;
        size_t _number_of_islands = 5;
        island_structure _island_structure = island_structure::random;
        island_migration_policy _island_migration_policy = island_migration_policy::random;
        island_replacement_policy _island_replacement_policy = island_replacement_policy::worst_swap;
        size_t _migration_epoch = 25;
        size_t _migration_size = 2;
        double _fitness_sharing_niche_size = 5.0;
        // Stopping criteria
        size_t _max_generations = 150;
        // Reproduction
        const size_t _parents_per_children = 2;
        double _children_proportion = 7.0;
        double _crossover_probability = 0.9;
        double _mutation_strength = 0.1;
        // Selection
        bool _competition_between_parents_and_children = false;
        unsigned int _tournament_size = 2;
        double _overselection_proportion = 0.1;
        size_t _roundrobin_tournament_size = 10;
        double _elitism_proportion = 0.01;
        // scaling
        double _sigma_bias = 1.0;
        double _sigma_constant = 2.0;
        double _linear_rank_selective_pressure = 1.5;

        // Data
        population_type _population;
        problem _problem;
        size_t _current_generation = 0;

        // Auxiliary generator
        static std::default_random_engine _generator;

        // Solution comparing
        std::function<bool(solution_ptr &, solution_ptr &)> _comp;
        std::function<bool(solution_ptr &, solution_ptr &)> _not_comp;
        std::function<bool(solution_ptr &, solution_ptr &)> _comp_fitness;
        std::function<bool(solution_ptr &, solution_ptr &)> _not_comp_fitness;

        // Get statistics
        population_type _best_solutions;
        size_t _number_of_best_solutions = 10;

};

#endif //EVOLUTIONARY_COMPUTATION_EA_H
