//
// Created by sam on 28/01/2022.
//

#include <libalgebra/libalgebra.h>
#include <libalgebra/operators.h>

#include <UnitTest++.h>
#include <libalgebra/alg_types.h>

SUITE(test_operator_smul) {

    struct fixture : alg_types<5, 5, Rational>
    {
        using multiplier_t = alg::operators::left_multiplication_operator<TENSOR>;

        using smul_operator_t = alg::operators::scalar_multiply_operator<multiplier_t, S>;
        using operator_t = alg::operators::linear_operator<smul_operator_t>;
        operator_t op;

        fixture() : op(multiplier_t(TENSOR(LET(1), S(2))), S(3))
        {}

    };

    TEST_FIXTURE(fixture, test_zero) {
        TENSOR zero, expected;

        CHECK_EQUAL(expected, op(zero));
    }

    TEST_FIXTURE(fixture, test_unit) {
        TENSOR unit(S(1)), expected;

        expected.add_scal_prod(typename TENSOR::KEY(LET(1)), S(6));

        CHECK_EQUAL(expected, op(unit));
    }

    TEST_FIXTURE(fixture, test_letter_1) {
        TENSOR arg(LET(1), S(1)), expected;

        expected.add_scal_prod(typename TENSOR::KEY{LET(1), LET(1)}, S(6));

        CHECK_EQUAL(expected, op(arg));
    }

    TEST_FIXTURE(fixture, test_letters_2345) {
        TENSOR arg, expected;

        arg.add_scal_prod(typename TENSOR::KEY(LET(2)), S(2));
        arg.add_scal_prod(typename TENSOR::KEY(LET(3)), S(3));
        arg.add_scal_prod(typename TENSOR::KEY(LET(4)), S(4));
        arg.add_scal_prod(typename TENSOR::KEY(LET(5)), S(5));

        expected.add_scal_prod(typename TENSOR::KEY {LET(1), LET(2)}, S(12));
        expected.add_scal_prod(typename TENSOR::KEY {LET(1), LET(3)}, S(18));
        expected.add_scal_prod(typename TENSOR::KEY {LET(1), LET(4)}, S(24));
        expected.add_scal_prod(typename TENSOR::KEY {LET(1), LET(5)}, S(30));
    }



}
