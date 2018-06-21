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
    class individual : public solution {
    public:
        individual(problem& rhs) : solution(rhs) {}
        individual(solution rhs) : solution(rhs) {}
        individual(solution& rhs) : solution(rhs) {}
        void mutation(problem &p, double mutation_strength);
        individual crossover(problem &p, individual& rhs);
        double _second_order_mutation_strength = 0.5;
        double _crossover_probability = 0.8;
        double _mutation_strength = 0.1;
    private:
        void reflect(double& v, double lb, double ub);
    };

    using solution_ptr = std::shared_ptr<solution>;
    using individual_ptr = std::shared_ptr<EA::individual>;
    using population_type = std::vector<individual_ptr>;

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

    enum class algorithm {
        DEFAULT, // Default configuration
        GA, // Genetic Algorithm
        EE // Evolutionary Strategy
    };

public:
    EA(problem &p);

    void algorithm(enum algorithm alg);

    void run();

    void run(size_t iterations);

    double best_fx();

    individual_ptr best_solution();

    typename population_type::iterator best_solutions_begin();

    typename population_type::iterator best_solutions_end();

    // Setters and getters
    // Population management
    size_t population_size();
    void population_size(size_t value);
    size_t number_of_islands();
    void number_of_islands(size_t value);
    island_structure island_structure();
    void island_structure(enum island_structure value);
    island_migration_policy island_migration_policy();
    void island_migration_policy(enum island_migration_policy value);
    island_replacement_policy island_replacement_policy();
    void island_replacement_policy(enum island_replacement_policy value);
    size_t migration_epoch();
    void migration_epoch(size_t value);
    size_t migration_size();
    void migration_size(size_t value);
    double fitness_sharing_niche_size();
    void fitness_sharing_niche_size(double value);
    // Stopping criteria
    size_t max_generations();
    void max_generations(size_t value);
    // Reproduction
    double children_proportion();
    void children_proportion(double value);
    double crossover_probability();
    void crossover_probability(double value);
    double mutation_strength();
    void mutation_strength(double value);
    // Selection
    bool competition_between_parents_and_children();
    void competition_between_parents_and_children(bool value);
    unsigned int tournament_size();
    void tournament_size(unsigned int value);
    double overselection_proportion();
    void overselection_proportion(double value);
    size_t roundrobin_tournament_size();
    void roundrobin_tournament_size(size_t value);
    double elitism_proportion();
    void elitism_proportion(double value);
    // scaling
    enum scaling_strategy reproduction_scaling_strategy();
    void reproduction_scaling_strategy(enum scaling_strategy value);
    enum selection_strategy reproduction_selection_strategy();
    void reproduction_selection_strategy(enum selection_strategy value);
    enum scaling_strategy survival_scaling_strategy();
    void survival_scaling_strategy(enum scaling_strategy value);
    enum selection_strategy survival_selection_strategy();
    void survival_selection_strategy(enum selection_strategy value);
    double sigma_bias();
    void sigma_bias(double value);
    double sigma_constant();
    void sigma_constant(double value);
    double linear_rank_selective_pressure();
    void linear_rank_selective_pressure(double value);

private /* methods */:
    // main cycle
    void evolutionary_cycle();

    // helpers
    size_t n_of_selection_candidates();
    size_t size_of_elite_set();
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
    population_type reproduction(population_type &parents);

    population_type merge(population_type &parents,population_type &children);

    // update
    population_type insert_elite_set(population_type& individuals, population_type& elite_source, size_t size_of_elite_set);

    void population_update_step(population_type &children);

    void initialize_population();

    void attribute_fitness_from_rank(population_type &population);

    population_type get_island(population_type &population, int idx);

    void try_to_update_best(individual_ptr &candidate);

    void migration_step();

private /* members */:
    // Search parameters
    // Population management
    size_t _population_size = 300;
    size_t _number_of_islands = 5;
    enum island_structure _island_structure = island_structure::random;
    enum island_migration_policy _island_migration_policy = island_migration_policy::random;
    enum island_replacement_policy _island_replacement_policy = island_replacement_policy::worst_swap;
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
    enum scaling_strategy _reproduction_scaling_strategy = scaling_strategy::window;
    enum selection_strategy _reproduction_selection_strategy = selection_strategy::tournament;
    enum scaling_strategy _survival_scaling_strategy = scaling_strategy::window;
    enum selection_strategy _survival_selection_strategy = selection_strategy::truncate;
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
    std::function<bool(individual_ptr &, individual_ptr &)> _comp;
    std::function<bool(individual_ptr &, individual_ptr &)> _not_comp;
    std::function<bool(individual_ptr &, individual_ptr &)> _comp_fitness;
    std::function<bool(individual_ptr &, individual_ptr &)> _not_comp_fitness;

    // Get statistics
    population_type _best_solutions;
    size_t _number_of_best_solutions = 10;

};

#endif //EVOLUTIONARY_COMPUTATION_EA_H
