#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <algorithm>
#include <type_traits>

#include "tsp.h"
#include "route.h"

#include "knapsack.h"
#include "knapsackP.h"

#include "real_function.h"
#include "real_vector.h"

template <typename problem_t, typename solution_t>
void random_search() {
    // Problem parameters
    const size_t problem_size = 20;

    // Create problem
    problem_t problem(problem_size);
    problem.disp();

    // Search parameters
    const size_t max_iterations = 100;

    // Get statistics
    double max = std::numeric_limits<double>::min();
    double min = std::numeric_limits<double>::max();
    double avg = 0.0;

    // Random search
    for (size_t i = 0; i < max_iterations; ++i) {
        solution_t candidate_solution(problem);
        std::cout << "Random solution #" << i + 1 << std::endl;
        candidate_solution.disp(problem);
        double f = candidate_solution.evaluate(problem);
        std::cout << "Fitness: " << f << std::endl << std::endl;
        avg += f;
        max = std::max(max,f);
        min = std::min(min,f);
    }

    // Print final statistics
    std::cout << "Max fx for " << typeid(problem).name() << ": " << max << std::endl;
    std::cout << "Avg fx for " << typeid(problem).name() << ": " << avg/max_iterations << std::endl;
    std::cout << "Min fx for " << typeid(problem).name() << ": " << min << std::endl;
}

int main() {
    random_search<tsp,route>();
    random_search<knapsack_p,knapsack>();
    random_search<real_function,real_vector>();
    return 0;
}