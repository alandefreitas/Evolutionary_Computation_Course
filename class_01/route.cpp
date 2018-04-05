//
// Created by Alan de Freitas on 05/04/2018.
//

#include "route.h"


std::default_random_engine route::_generator = std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count());

route::route(tsp &p)
        : _route(p.size())
{
    std::iota(_route.begin(), _route.end(), 0);
    std::shuffle(_route.begin(), _route.end(), this->_generator);
}

void route::disp(tsp &p) {
    std::cout << "Solution" << std::endl;
    for (int i = 0; i < p.size(); ++i) {
        std::cout << this->_route[i] << "\t";
    }
    std::cout << std::endl;
}

double route::evaluate(tsp &p) {
    double total = 0.0;
    for (int i = 0; i < p.size() - 1; ++i) {
        total += p.distance(this->_route[i],this->_route[i+1]);
    }
    total += p.distance(this->_route[p.size() - 1],this->_route[0]);
    return total;
}