//
// Created by Alan de Freitas on 05/04/2018.
//

#include "tsp.h"

std::default_random_engine tsp::_generator = std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count());

tsp::tsp(size_t n)
        : _distance(n,std::vector<double>(n,0.0))
{
    std::uniform_real_distribution<double> d(100.00, 500.00);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            this->_distance[i][j] = d(this->_generator);
            this->_distance[j][i] = this->_distance[i][j];
        }
    }
}

void tsp::disp() {
    std::cout << "Problem" << std::endl;
    for (int i = 0; i < this->size(); ++i) {
        for (int j = 0; j < this->size(); ++j) {
            std::cout << this->_distance[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}

size_t tsp::size() {
    return this->_distance.size();
}

double tsp::distance(size_t i, size_t j) {
    return this->_distance[i][j];
}
