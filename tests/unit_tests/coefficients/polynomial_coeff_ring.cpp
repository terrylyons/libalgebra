//
// Created by sam on 18/07/22.
//

#include "libalgebra/polynomial_coefficients.h"
#include <UnitTest++.h>
#include <libalgebra/implementation_types.h>

using namespace alg;

SUITE(polynomial_ring_tests)
{

    struct fixture {
        static constexpr DEG WIDTH = 5;
        static constexpr DEG DEPTH = 5;

        using coeff_ring = coefficients::rational_poly_ring;
        using rational = typename coeff_ring::RAT;

        using S = typename coeff_ring::SCA;
        using Q = typename coeff_ring::RAT;
        using key_type = typename S::KEY;
    };

    TEST_FIXTURE(fixture, test_addition)
    {
        S one(LET(1), rational(1));
        S two(LET(2), rational(1));

        S expected = one + two;

        CHECK_EQUAL(expected, coeff_ring::add(one, two));
    }
}
