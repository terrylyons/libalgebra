//
// Created by sam on 24/01/2022.
//
#include "standard_coefficient_test_suite.h"
#include "unital_coefficient_test_suite.h"
#include <libalgebra/coefficients.h>

MAKE_SUITE_FOR(float_field_standard_test_suite, standard_coefficient_tests, alg::coefficients::float_field)
MAKE_SUITE_FOR(float_field_unital_test_suite, unital_coefficients_tests, alg::coefficients::float_field)
