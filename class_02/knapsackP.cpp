//
// Created by Alan de Freitas on 05/04/2018.
//

#include "knapsackP.h"

std::default_random_engine knapsack_p::_generator = std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count());

knapsack_p::knapsack_p(size_t n)
        : _value(n,0),
          _weight(n,0),
          _capacity(0)
{
    std::uniform_int_distribution<int> d(100, 500);
    for (int i = 0; i < n; ++i) {
        this->_value[i] = d(this->_generator);
        this->_weight[i] = d(this->_generator);
    }
    this->_capacity = std::accumulate(this->_weight.begin(),this->_weight.end(),0)/2;
    this->_max_value = std::accumulate(this->_value.begin(),this->_value.end(),0);
}

void knapsack_p::disp() {
    std::cout << "Knapsack Problem" << std::endl;
    std::cout << "Value: ";
    for (int i = 0; i < this->size(); ++i) {
        std::cout << this->_value[i];
    }
    std::cout << std::endl;
    std::cout << "Weight: ";
    for (int i = 0; i < this->size(); ++i) {
        std::cout << this->_weight[i];
    }
    std::cout << std::endl;
    std::cout << "Capacity: " << this->_capacity << std::endl;
}

size_t knapsack_p::size() {
    return this->_value.size();
}

int knapsack_p::value(size_t i){
    return this->_value[i];
}

int knapsack_p::weight(size_t i){
    return this->_weight[i];
}

int knapsack_p::capacity(){
    return this->_capacity;
}

int knapsack_p::max_value(){
    return this->_max_value;
}