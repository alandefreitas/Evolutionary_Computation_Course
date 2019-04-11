#include <iostream>

// knapsack
#include "knapsackP.h"
#include "knapsack.h"
// TSP
#include "route.h"
#include "tsp.h"
// real function
#include "real_function.h"
#include "real_vector.h"
// EA
#include "naive_EA.hpp"

template <typename problem_class, typename solution_class>
void solve_problem(size_t problem_size){
    // Create a random problem
    problem_class problem(problem_size);
    problem.disp();

    // Create a solver for the problem
    naive_EA<problem_class,solution_class> solver(problem);
    solver.run();

    // Print final statistics
    std::cout << "Best solution: " << solver.best_fx() << std::endl;
};

int main() {
    const size_t problem_size = 20;
    solve_problem<knapsack_p,knapsack>(problem_size);
    solve_problem<real_function,real_vector>(problem_size);
    solve_problem<tsp,route>(problem_size);
    return 0;
}