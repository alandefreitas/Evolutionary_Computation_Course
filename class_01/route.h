//
// Created by Alan de Freitas on 05/04/2018.
//

#ifndef EVOLUTIONARY_COMPUTATION_ROUTE_H
#define EVOLUTIONARY_COMPUTATION_ROUTE_H

#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <algorithm>

#include "tsp.h"

class route {
    public:
        // Functions we need for a random search
        explicit route(tsp &p);

        void disp(tsp &p);

        double evaluate(tsp &p);

        // Functions we need for an evolutionary algorithm
        void mutation(tsp &p, double mutation_strength);

        route crossover(tsp &p, route &rhs);

        // Auxiliary functions
        double distance(tsp &p, route &rhs, double max_dist = std::numeric_limits<double>::max());

    private:
        // Private members
        std::vector<size_t> _route;
        static std::default_random_engine _generator;
};


#endif //EVOLUTIONARY_COMPUTATION_ROUTE_H
