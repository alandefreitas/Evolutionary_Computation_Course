//
// Created by Alan de Freitas on 05/04/2018.
//

#include "real_vector.h"


std::default_random_engine real_vector::_generator = std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count());

real_vector::real_vector(real_function &p)
        : _real_vector(p.size())
{
    std::uniform_real_distribution<double> d(p.lower_bound(),p.upper_bound());
    for (double& item : this->_real_vector) {
        item = d(this->_generator);
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
    const double A = 10;
    const double n = p.size();
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += pow(_real_vector[i],2) - A * cos(2. * 3.14 * _real_vector[i]);
    }
    return A * n + sum;
}

void real_vector::mutation(real_function &p, double mutation_strength) {
    std::normal_distribution<double> d(0,mutation_strength);
    for (double& item : this->_real_vector) {
        item += d(_generator);
        while (item < p.lower_bound() || item > p.upper_bound()) {
            if (item < p.lower_bound()) {
                item = item + 2 * (p.lower_bound() - item);
            } else if (item > p.upper_bound()) {
                item = item - 2 * (item - p.upper_bound());
            }
        }
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
