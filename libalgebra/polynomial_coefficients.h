//
// Created by sam on 18/07/22.
//

#ifndef LIBALGEBRA_POLYNOMIAL_COEFFICIENTS_H
#define LIBALGEBRA_POLYNOMIAL_COEFFICIENTS_H

#include "coefficients.h"
#include "polynomials.h"
#include "rational_coefficients.h"

namespace alg {
namespace coefficients {

template<typename BaseCoefficientRing>
using poly_ring = coefficient_ring<poly<BaseCoefficientRing>, typename BaseCoefficientRing::RAT>;

using float_poly_ring = poly_ring<float_field>;
using double_poly_ring = poly_ring<double_field>;
using rational_poly_ring = poly_ring<rational_field>;

}// namespace coefficients
}// namespace alg

#endif//LIBALGEBRA_POLYNOMIAL_COEFFICIENTS_H
