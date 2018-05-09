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

void route::mutation(tsp &p, double mutation_strength) {
    std::uniform_int_distribution<size_t> d(0,p.size()-1);
    for (int i = 0; i < (mutation_strength / 2) * p.size(); ++i) {
        std::swap(this->_route[d(_generator)],this->_route[d(_generator)]);
    }
}

route route::crossover(tsp &p, route &rhs) {
    std::uniform_int_distribution<size_t> d(0,p.size()-1);
    route child(p);
    std::vector<int> set(p.size(),0);
    // copy from parent 1
    size_t pos1 = d(_generator);
    size_t pos2 = d(_generator);
    if (pos1 > pos2) {
        std::swap(pos1,pos2);
    }
    std::copy(this->_route.begin()+pos1,
              this->_route.begin()+pos2,
              child._route.begin()+pos1);
    for (int i = pos1; i < pos2; ++i) {
        set[this->_route[i]] = 1;
    }
    // copy from parent2
    size_t k = pos2;
    for (int i = pos2; i < p.size(); ++i) {
        if (!set[rhs._route[i]]){
            child._route[k % p.size()] = rhs._route[i];
            k++;
        }
    }
    for (int i = 0; i < pos2; ++i) {
        if (!set[rhs._route[i]]){
            child._route[k % p.size()] = rhs._route[i];
            k++;
        }
    }
    return child;
}

double route::distance(tsp &p, route &rhs, double max_dist) {
    double hamming = 0.0;
    for (int i = 0; i < p.size(); ++i) {
        const size_t this_city = this->_route[i%p.size()];
        size_t next_city = this->_route[(i+1)%p.size()];
        auto iter = std::find(rhs._route.begin()+i,rhs._route.end(),this_city);
        if (iter == rhs._route.end()){
            iter = std::find(rhs._route.begin(),rhs._route.end(),this_city);
        }
        size_t this_city_pos = iter - rhs._route.begin();
        size_t next_city_rhs = rhs._route[(this_city_pos+1)%p.size()];
        if (next_city != next_city_rhs){
            hamming++;
            if (hamming > max_dist){
                return hamming;
            }
        }
    }
    return hamming;
}
