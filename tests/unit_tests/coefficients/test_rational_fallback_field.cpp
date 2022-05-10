//
// Created by sam on 24/01/2022.
//
#include "libalgebra/coefficients/coefficients.h"

#include <boost/multiprecision/cpp_int.hpp>

#include "standard_coefficient_test_suite.h"
#include "unital_coefficient_test_suite.h"

namespace mp = boost::multiprecision;

using number = mp::number<mp::cpp_rational_backend, mp::et_off>;
using reserve_field = alg::coefficients::coefficient_field<number>;

MAKE_SUITE_FOR(rational_fallback_field_standard_test_suite, standard_coefficient_tests, reserve_field);
MAKE_SUITE_FOR(rational_fallback_field_unital_test_suite, unital_coefficients_tests, reserve_field);
