#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <algorithm>
#include <type_traits>

#include "tsp.h"
#include "route.h"

int main() {
    // Problem parameters
    const size_t problem_size = 20;

    // Create problem
    tsp problem(problem_size);
    problem.disp();

    // Search parameters
    const size_t max_iterations = 10000000;

    // Get statistics
    double max = std::numeric_limits<double>::min();
    double min = std::numeric_limits<double>::max();
    double avg = 0.0;

    // Random search
    for (size_t i = 0; i < max_iterations; ++i) {
        route candidate_solution(problem);
        std::cout << "Random solution #" << i + 1 << std::endl;
        candidate_solution.disp(problem);
        double f = candidate_solution.evaluate(problem);
        std::cout << "Fitness: " << f << std::endl << std::endl;
        avg += f;
        max = std::max(max,f);
        min = std::min(min,f);
    }

    // Print final statistics
    std::cout << "Max route size: " << max << std::endl;
    std::cout << "Avg route size: " << avg/max_iterations << std::endl;
    std::cout << "Min route size: " << min << std::endl;

    return 0;
}