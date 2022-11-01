//
// Created by sam on 02/02/2021.
//

#include <UnitTest++.h>

#include <libalgebra/libalgebra.h>
#include <libalgebra/alg_types.h>

#include "../../common/time_and_details.h"
#include "../../common/helpers.h"


using alg::DEG;
using alg::LET;

struct Fixture : alg_types<5, 5, DPReal>
{
    typedef alg_types<5, 5, DPReal> AT;
    typedef typename AT::S S;
    typedef typename AT::Q Q;
    typedef typename AT::TENSOR TENSOR;
    typedef typename TENSOR::KEY TKEY;
};



SUITE(test_tensor_functions) {

    template<class T, class S>
    T exp_to_depth(T x, size_t depth, S one) {
        T result (one);
        T xn(one);
        S fact = 1;
        for (size_t i=1; i<=depth; ++i) {
            xn *= x;
            fact *= double(i);
            result += (xn / fact);
        }
        return result;
    }

    TEST_FIXTURE(Fixture, tensor_exponential_zero) {
        TEST_DETAILS();
        TENSOR ten(S(0));
        TENSOR expected (S(1));

        CHECK_EQUAL(expected, exp(ten));
    }

    TEST_FIXTURE(Fixture, test_exponential_tensor_unit) {
        TEST_DETAILS();

        TENSOR ten (TENSOR::VECT::one);
        TENSOR expected(exp_to_depth(S(1), 5, S(1)));
        TENSOR result(exp(ten));
        CHECK_VEC_CLOSE(expected, result, 2.0e-15);
    }

    TEST_FIXTURE(Fixture, test_exponential_single_letter) {
        TEST_DETAILS();
        LET letter = 1;

        TENSOR ten (letter, S(1));
        TENSOR expected = exp_to_depth(ten, 5, S(1));
        TENSOR result (exp(ten));

        CHECK_VEC_CLOSE(expected, result, 2.0e-15);
    }

    TEST_FIXTURE(Fixture, test_exponential_multiple_letter) {
        TEST_DETAILS();

        TENSOR ten (LET(1), S(1));
        ten += TENSOR(LET(2), S(1));
        TENSOR expected = exp_to_depth(ten, 5, S(1));

        CHECK_VEC_CLOSE(expected, exp(ten), 2.0e-15);
    }

    TEST_FIXTURE(Fixture, test_exponential_single_letter_with_coeff) {
        TEST_DETAILS();
        LET letter = 1;

        TENSOR ten (letter, S(2));
        TENSOR expected = exp_to_depth(ten, 5, S(1));

        CHECK_VEC_CLOSE(expected, exp(ten), 2.0e-15);
    }

    TEST_FIXTURE(Fixture, test_exponential_multiple_letter_with_coeffs) {
        TEST_DETAILS();

        TENSOR ten (LET(1), S(1));
        ten += TENSOR(LET(2), S(2));
        TENSOR expected = exp_to_depth(ten, 5, S(1));

        CHECK_VEC_CLOSE(expected, exp(ten), 2.0e-15);
    }

    TEST_FIXTURE(Fixture, test_log_tensor_unit) {
        TEST_DETAILS();

        TENSOR tunit ( TENSOR::VECT::one );
        TENSOR zero (S(0));

        CHECK_EQUAL(zero, log(tunit));
    }

    TEST_FIXTURE(Fixture, test_log_tensor_with_no_explicit_unit) {
        TEST_DETAILS();

        TENSOR tunit ( TENSOR::VECT::one );
        TENSOR no_unit (LET(1), S(1));
        TENSOR ten = tunit + no_unit;

        CHECK_EQUAL(log(ten), log(no_unit));
    }

    TEST_FIXTURE(Fixture, test_log_exp_round_trip_single_letter) {
        TEST_DETAILS();

        TENSOR ten (LET(1), S(1));

        CHECK_VEC_CLOSE(ten, log(exp(ten)), 2.0e-15);
    }

    TEST_FIXTURE(Fixture, test_exp_log_round_trip_single_letter) {
        TEST_DETAILS();

        TENSOR ten (LET(1), S(1));
        TENSOR expected = TENSOR(TKEY()) + ten;

        CHECK_VEC_CLOSE(expected, exp(log(ten)), 2.0e-15);
    }




}
