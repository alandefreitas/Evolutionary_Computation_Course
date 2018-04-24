//
// Created by Alan de Freitas on 05/04/2018.
//

#include "knapsack.h"


std::default_random_engine knapsack::_generator = std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count());

knapsack::knapsack(knapsack_p &p)
        : _knapsack(p.size())
{
    std::uniform_int_distribution<int> d(0,1);
    for (int& item : this->_knapsack) {
        item = d(this->_generator);
    }
}

void knapsack::disp(knapsack_p &p) {
    std::cout << "Solution" << std::endl;
    for (int i = 0; i < p.size(); ++i) {
        std::cout << this->_knapsack[i] << "\t";
    }
    std::cout << std::endl;
}

double knapsack::evaluate(knapsack_p &p) {
    int total_value = 0.0;
    int total_weight = 0.0;
    for (int i = 0; i < p.size() - 1; ++i) {
        total_value += this->_knapsack[i] * p.value(i);
        total_weight += this->_knapsack[i] * p.weight(i);
    }
    if (total_weight > p.capacity()){
        total_value -= p.max_value() * (total_weight - p.capacity());
    }
    return total_value;
}

void knapsack::mutation(knapsack_p &p, double mutation_strength) {
    std::uniform_real_distribution<double> d(0.0,1.0);
    for (int& item : this->_knapsack) {
        if (d(_generator) < mutation_strength) {
            item = 1 - item;
        }
    }
    // Smart approach:
    // std::binomial_distribution<double> d(p.size(), mutation_strength)
    // int n = d(_generator)
    // flip $n$ bits
}

knapsack knapsack::crossover(knapsack_p &p, knapsack &rhs) {
    std::uniform_int_distribution<int> d(0,1);
    knapsack child(p);
    for (int i = 0; i < p.size(); ++i) {
        if (d(_generator)) {
            child._knapsack[i] = this->_knapsack[i];
        } else {
            child._knapsack[i] = rhs._knapsack[i];
        }
    }
    return child;
}
