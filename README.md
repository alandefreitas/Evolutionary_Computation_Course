# Evolutionary Computation (Examples)

Examples for the course on evolutionary computation:

* [Class 01](./class_01/): Modelling and Random Search 
* [Class 02](./class_02/): Representation and Variation Operators
    * Binary
    * Integer
    * Real
    * Permutations
    * Tree
* [Class 03](./class_03/): Selection and Scaling
    * Selection
        * Uniform (same probability to all solutions)
        * Truncate (only the solutions with best fitness)
        * Tournament (winners of tournaments between groups of solutions)
        * Roulette (selection probability proportional to fitness)
        * Stochastic Uniform Selection (stochastic probability proportional to fitness)
        * Overselection (80% of operations on the x% best solutions)
        * Round-Robin Tournament (each individual evaluated against other q individuals)
    * Scaling
        * Window (`fitness <- fx - worst_fx`)
        * Sigma (`fitness <- max(bias + (fx - avg_fx)/(c * std_fx)), 0.0)`)
        * Linear rank (`fitness <- ((2-pressure_constant)/n_ind) + ((2 * rank(fx) * (pressure_constant - 1))/n_ind * (n_ind - 1))`)
        * Exponential rank (`fitness <- (1-e^(-i))/n_ind`)

## Running the examples


I have been posting some examples from my course on Evolutionary Computation here.

They are examples on how to define a new problem and solve it with an Evolutionary Algorithm.
 
There’s a basic EA and I've been posting new examples of more components every 1 or 2 weeks. 

All you have to do is create new classes using `KnapsackP`/`Knapsack` or `TSP`/`route` as reference. It'd be nice if you have new problems to contribute.

The code uses very basic C++. It should work in any compiler with support to *C++11*. Just put the files in a new C++ project and you’re good to go.

The project is organized with CMake, which you can download from [their website](www.cmake.org). Don’t download from repositories because they are usually very old versions of CMake. 

Once you have CMake installed all you have to do is go to the terminal run `cmake` on the project folder.

If you have an IDE that supports CMake (such as CLion), it’s even easier. The process should be transparent. Just push *play*.

However, these examples are didactic. We are now working on a robust framework/solver based on Evolutionary Algorithms and Symbolic Computation, which is needed for difficult large scale dynamic problems. The source code, however, is not available yet because we are still writing the white paper.