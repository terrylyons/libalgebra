//
// Created by sam on 24/01/2022.
//
#include "standard_coefficient_test_suite.h"
#include "unital_coefficient_test_suite.h"
#include <libalgebra/coefficients.h>
#include <libalgebra/mpfloat_coefficients.h>

MAKE_SUITE_FOR(mpfloat_field_standard_test_suite, standard_coefficient_tests, alg::coefficients::mpfloat_field<20U>)
MAKE_SUITE_FOR(mpfloat_field_unital_test_suite, unital_coefficients_tests, alg::coefficients::mpfloat_field<20U>)
