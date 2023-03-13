//
// Created by sam on 24/01/2022.
//

#include "standard_coefficient_test_suite.h"
#include "unital_coefficient_test_suite.h"
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <libalgebra/coefficients.h>

namespace mp = boost::multiprecision;

template<DEG Digits10>
using number = mp::number<mp::cpp_bin_float<Digits10>, mp::et_off>;

template<DEG Digits10>
using reserve_field = alg::coefficients::coefficient_field<number<Digits10>>;

MAKE_SUITE_FOR(float_mpfield_fallback_standard_test_suite, standard_coefficient_tests, reserve_field<20U>)
MAKE_SUITE_FOR(float_mpfield_fallback_unital_test_suite, unital_coefficients_tests, reserve_field<20U>)
