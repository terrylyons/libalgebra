//
// Created by sam on 28/01/2022.
//

#include <libalgebra/libalgebra.h>
#include <libalgebra/operators.h>

#include <UnitTest++.h>
#include <libalgebra/alg_types.h>

SUITE(test_operator_sum)
{

    struct fixture : alg_types<5, 5, Rational> {

        using multiplier_t = alg::operators::left_multiplication_operator<TENSOR>;

        using sum_type = alg::operators::sum_operator<multiplier_t, multiplier_t>;

        using operator_t = alg::operators::linear_operator<sum_type>;

        static multiplier_t make_left_tensor()
        {
            return multiplier_t(TENSOR(LET(1), S(1)));
        }

        static multiplier_t make_right_tensor()
        {
            return multiplier_t(TENSOR(LET(2), S(2)));
        }

        operator_t op;

        fixture() : op(make_left_tensor(), make_right_tensor())
        {}
    };

    TEST_FIXTURE(fixture, test_zero)
    {
        TENSOR zero, expected;

        CHECK_EQUAL(expected, op(zero));
    }

    TEST_FIXTURE(fixture, test_unit)
    {
        TENSOR unit(S(1));

        TENSOR expected;
        expected.add_scal_prod(typename TENSOR::KEY(LET(1)), S(1));
        expected.add_scal_prod(typename TENSOR::KEY(LET(2)), S(2));

        CHECK_EQUAL(expected, op(unit));
    }

    TEST_FIXTURE(fixture, test_letter_1) {
        TENSOR arg(LET(1), S(1));

        TENSOR expected;
        expected.add_scal_prod(typename TENSOR::KEY{LET(1), LET(1)}, S(1));
        expected.add_scal_prod(typename TENSOR::KEY{LET(2), LET(1)}, S(2));

        CHECK_EQUAL(expected, op(arg));
    }
}
