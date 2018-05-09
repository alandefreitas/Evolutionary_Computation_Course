#include <iostream>

#include "knapsackP.h"
#include "knapsack.h"
#include "route.h"
#include "tsp.h"
#include "EA.hpp"

template <typename problem_class, typename solution_class>
void solve_problem(size_t problem_size){
    // Create a random problem
    problem_class problem(100);
    problem.disp();

    // Create a solver for the problem
    EA<problem_class,solution_class> solver(problem);
    solver.run();

    // Print final statistics
    int idx = 1;
    for (auto iter = solver.best_solutions_begin(); iter != solver.best_solutions_end(); ++iter) {
        std::cout << "Solution " << idx << ": " << std::endl;
        (*iter)->disp(problem);
        std::cout << "Objetive function " << idx << ": " << (*iter)->fx << std::endl;
        idx++;
    }
};

int main() {
    const size_t problem_size = 20;
    solve_problem<knapsack_p,knapsack>(problem_size);
    solve_problem<tsp,route>(problem_size);
    return 0;
}