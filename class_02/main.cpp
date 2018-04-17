#include <iostream>

#include "knapsackP.h"
#include "knapsack.h"
#include "route.h"
#include "tsp.h"
#include "EA.hpp"

template <typename problem_class, typename solution_class>
void solve_problem(size_t problem_size){
    // Create a random problem
    problem_class problem(20);
    problem.disp();

    // Create a solver for the problem
    EA<problem_class,solution_class> s(problem);
    s.run();

    // Print final statistics
    std::cout << "Best solution: " << s.best_fx() << std::endl;
};

int main() {
    const size_t problem_size = 20;
    solve_problem<knapsack_p,knapsack>(problem_size);
    solve_problem<tsp,route>(problem_size);
    return 0;
}