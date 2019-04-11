//
// Created by Alan de Freitas on 05/04/2018.
//

#include "real_vector.h"


std::default_random_engine real_vector::_generator = std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count());

real_vector::real_vector(real_function &p)
        : _real_vector(p.size())
{
    for (size_t i = 0; i < this->_real_vector.size(); ++i) {
        std::uniform_real_distribution<double> d(p.lower_bound(i),p.upper_bound(i));
        _real_vector[i] = d(this->_generator);
    }
}

void real_vector::disp(real_function &p) {
    std::cout << "Solution" << std::endl;
    for (int i = 0; i < p.size(); ++i) {
        std::cout << this->_real_vector[i] << "\t";
    }
    std::cout << std::endl;
}

double real_vector::evaluate(real_function &p) {
    return p.fx(this->_real_vector);
}

void reflection(double& value, double lower_bound, double upper_bound) {
    while (value < lower_bound || value > upper_bound) {
        if (value < lower_bound) {
            value = value + 2 * (lower_bound - value);
        } else if (value > upper_bound) {
            value = value - 2 * (value - upper_bound);
        }
    }
}

void real_vector::mutation(real_function &p, double mutation_strength) {
    std::normal_distribution<double> d(0,mutation_strength);
    for (size_t i = 0; i < this->_real_vector.size(); ++i) {
        this->_real_vector[i] += d(_generator);
        reflection(this->_real_vector[i], p.lower_bound(i), p.upper_bound(i));
    }
}

real_vector real_vector::crossover(real_function &p, real_vector &rhs) {
    std::uniform_real_distribution<double> d(0.,1.);
    double alpha = d(_generator);
    real_vector child(p);
    for (int i = 0; i < p.size(); ++i) {
        child._real_vector[i] = alpha * this->_real_vector[i] + (1-alpha) * rhs._real_vector[i];
    }
    return child;
}

double real_vector::distance(real_function &p, real_vector &rhs, double max_dist) {
    double total = 0.0;
    for (int i = 0; i < rhs._real_vector.size(); ++i) {
        total += std::abs(rhs._real_vector[i] - _real_vector[i]);
        if (total > max_dist) {
            return total;
        }
    }
    return total;
}
