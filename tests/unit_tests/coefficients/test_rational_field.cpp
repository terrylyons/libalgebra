//
// Created by sam on 24/01/2022.
//

#include <boost/multiprecision/gmp.hpp>

#include "standard_coefficient_test_suite.h"
#include "unital_coefficient_test_suite.h"
#include <libalgebra/coefficients/coefficients.h>
#include <libalgebra/coefficients/rational_coefficients.h>

MAKE_SUITE_FOR(rational_field_standard_test_suite, standard_coefficient_tests, alg::coefficients::rational_field);
MAKE_SUITE_FOR(rational_field_unital_test_suite, unital_coefficients_tests, alg::coefficients::rational_field);
