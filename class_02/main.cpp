#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <algorithm>
#include <type_traits>

#include "knapsackP.h"
#include "knapsack.h"
#include "EA.hpp"

int main() {
    // Create problem
    knapsack_p problem(200);
    problem.disp();

    // Create solver for the knapsack problem
    EA<knapsack_p,knapsack> s(problem);
    s.run();

    // Print final statistics
    std::cout << "Max knapsack quality: " << s.max_fx() << std::endl;

    return 0;
}