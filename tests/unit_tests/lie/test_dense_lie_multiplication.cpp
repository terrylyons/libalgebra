//
// Created by sam on 05/02/2021.
//

/*
 * Currently, these tests don't work because the Lie basis mutliplication
 * is hard coded to output a LIE object with the default vector type rather
 * than a modified class.
 */
#if 0

#include <UnitTest++.h>

#include <libalgebra/libalgebra.h>
#include <libalgebra/alg_types.h>

#include "../time_and_details.h"


SUITE(dense_lie_multiplication) {


    struct Fixture : public alg_types<5, 5, Rational>
    {
        typedef alg_types<5, 5, Rational> ALG_TYPES;
        typedef typename ALG_TYPES::LIE::BASIS LBASIS;
        typedef typename LBASIS::KEY LKEY;

        struct field
        {
            typedef typename ALG_TYPES::S S;
            typedef typename ALG_TYPES::Q Q;
        };

        typedef alg::vectors::dense_vector<LBASIS, field> DENSE;
        typedef alg::algebra<LBASIS, field, DENSE> LIE;


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


#endif
