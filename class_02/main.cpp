#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <algorithm>
#include <type_traits>

#include "knapsackP.h"
#include "knapsack.h"
#include "EA.h"

int main() {
    // Problem parameters
    const size_t problem_size = 200;

    // Create problem
    knapsack_p problem(problem_size);
    problem.disp();

    // Create solver
    EA s(problem);
    s.run();

    // Print final statistics
    std::cout << "Max knapsack quality: " << s.max_fx() << std::endl;

    return 0;
}