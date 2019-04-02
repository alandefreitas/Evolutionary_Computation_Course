//
// Created by Alan de Freitas on 05/04/2018.
//

#include "real_function.h"


real_function::real_function(size_t n)
        : _size(n)
{
}

void real_function::disp() {
    std::cout << "Real function of size " << _size << std::endl;
}

size_t real_function::size() {
    return this->_size;
}