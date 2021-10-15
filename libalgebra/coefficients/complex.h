//
// Created by sam on 15/10/2021.
//

#ifndef LIBALGEBRA_COMPLEX_H
#define LIBALGEBRA_COMPLEX_H

#include <complex>

#include "coefficients.h"

namespace alg {
namespace coefficients {

template <typename Real=double>
using complex_field = coefficient_field<std::complex<Real>, std::complex<Real> >;


} // namespace coefficients
} // namsepace alg





#endif //LIBALGEBRA_COMPLEX_H
