//
// Created by sam on 05/02/2021.
//

#include <UnitTest++.h>

#include <libalgebra/libalgebra.h>
#include <libalgebra/alg_types.h>

#include "../../common/time_and_details.h"


SUITE(lie_multiplication) {


    struct Fixture : public alg_types<5, 5, Rational>
    {
        typedef alg_types<5, 5, Rational> ALG_TYPES;
        typedef typename ALG_TYPES::LIE LIE;
        typedef typename LIE::BASIS LBASIS;
        typedef typename LIE::KEY LKEY;

        const LIE lzero;
    };


    TEST_FIXTURE(Fixture, test_unidim_letter_by_zero) {
        TEST_DETAILS();

        LIE lhs(LKEY(1)), rhs(lzero);
        LIE expected(lzero);

        CHECK_EQUAL(expected, lhs * rhs);
    }

    TEST_FIXTURE(Fixture, test_unidim_letter_by_unidim_letter_ordered) {
        TEST_DETAILS();

        LIE lhs(LKEY(1)), rhs(LKEY(2));
        LIE expected(LKEY(6)); // 6 = [1,2]

        CHECK_EQUAL(expected, lhs * rhs);
    }

    TEST_FIXTURE(Fixture, test_unidim_letter_by_unidim_letter_unordered) {
        TEST_DETAILS();

        LIE lhs(LKEY(2)), rhs(LKEY(1));
        LIE expected = -LIE(LKEY(6)); // 6 = [1,2], so [2, 1] = -6

        CHECK_EQUAL(expected, lhs * rhs);
    }

    TEST_FIXTURE(Fixture, test_unidim_letter_by_unidim_non_letter_ordered) {
        TEST_DETAILS();

        LIE lhs(LKEY(1)), rhs(LKEY(6)); // 6 = [1,2]
        LIE expected(LKEY(16)); // 16 = [1,[1,2]]

        CHECK_EQUAL(expected, lhs * rhs);
    }

    TEST_FIXTURE(Fixture, test_unidim_product_degree_overflow) {
        TEST_DETAILS();

        LIE lhs(LKEY(16)), rhs(LKEY(17)); // 16, 17 both deg 3
        LIE expected(lzero);

        CHECK_EQUAL(expected, lhs * rhs);
    }

}
