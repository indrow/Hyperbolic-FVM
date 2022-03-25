//
// Created by indro on 10/03/21.
//

#ifndef HERTZ_CPI_H
#define HERTZ_CPI_H

#include <cmath>

template<typename real>
constexpr real CPi() {
    return std::atan(1.0) * 4.0;
}

#endif //HERTZ_CPI_H
