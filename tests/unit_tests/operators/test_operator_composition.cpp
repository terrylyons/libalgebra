//
// Created by sam on 28/01/2022.
//

#include <libalgebra/libalgebra.h>
#include <libalgebra/functionals.h>
#include <libalgebra/operators.h>

#include <UnitTest++.h>
#include <libalgebra/alg_types.h>

SUITE(operator_composition)
{

    struct fixture : alg_types<5, 5, Rational> {
        using types = alg_types<5, 5, Rational>;

        using multiply_t = alg::operators::left_multiplication_operator<TENSOR>;
        using functional_t = alg::operators::shuffle_tensor_functional<SHUFFLE_TENSOR, TENSOR>;

        using operator_t = alg::operators::linear_operator<
                alg::operators::composition_operator<functional_t, multiply_t>>;

        operator_t op;

        static functional_t make_functional()
        {
            SHUFFLE_TENSOR shuf;

            for (alg::LET l = 1; l <= ALPHABET_SIZE; ++l) {
                shuf.add_scal_prod(typename SHUFFLE_TENSOR::KEY(l), COEFF::one);
            }

            for (alg::LET l1 = 1; l1 <= ALPHABET_SIZE; ++l1) {
                for (alg::LET l2 = 1; l2 <= ALPHABET_SIZE; ++l2) {
                    shuf.add_scal_prod(typename SHUFFLE_TENSOR::KEY{l1, l2}, COEFF::one / S(2));
                }
            }

            return functional_t(shuf);
        }

        static multiply_t make_multiply()
        {
            TENSOR tens(alg::LET(1), S(2));
            return multiply_t(tens);
        }

        fixture() : op(make_multiply(), make_functional())
        {
        }
    };


    TEST_FIXTURE(fixture, test_zero) {
        TENSOR zero;

        CHECK_EQUAL(S(0), op(zero));
    }

    TEST_FIXTURE(fixture, test_unit_image) {
        TENSOR unit(S(1));

        CHECK_EQUAL(S(2), op(unit));
    }

    TEST_FIXTURE(fixture, test_image_letter_1) {
        TENSOR arg(alg::LET(1), S(1));

        CHECK_EQUAL(S(1), op(arg));
    }

    TEST_FIXTURE(fixture, test_image_letters_2345) {
        TENSOR arg;
        arg.add_scal_prod(typename TENSOR::KEY(alg::LET(2)), S(2));
        arg.add_scal_prod(typename TENSOR::KEY(alg::LET(3)), S(3));
        arg.add_scal_prod(typename TENSOR::KEY(alg::LET(4)), S(4));
        arg.add_scal_prod(typename TENSOR::KEY(alg::LET(5)), S(5));

        S expected = S(2) + S(3) + S(4) + S(5);

        CHECK_EQUAL(expected, op(arg));
    }

    TEST_FIXTURE(fixture, test_image_key11) {
        TENSOR arg(typename TENSOR::KEY { 1, 1});

        CHECK_EQUAL(S(0), op(arg));
    }


}
