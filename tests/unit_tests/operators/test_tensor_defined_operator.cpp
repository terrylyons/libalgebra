//
// Created by sam on 27/01/2022.
//

#include <libalgebra/alg_types.h>
#include <libalgebra/libalgebra.h>
#include <libalgebra/functionals.h>
#include <libalgebra/operators.h>
#include <libalgebra/tensor_operator.h>

#include <UnitTest++.h>

SUITE(tensor_defined_operator)
{

    struct fixture : alg_types<5, 5, Rational> {
        using types = alg_types<5, 5, Rational>;

        using out_types = alg_types<5, 4, Rational>;
        using OUT_TENSOR = typename out_types::TENSOR;
        using TKEY = typename TENSOR::KEY;
        using OKEY = typename OUT_TENSOR::KEY;

        using functional_t = alg::operators::shuffle_tensor_functional<SHUFFLE_TENSOR, TENSOR>;
        using op_impl = alg::operators::tensor_defined_operator<
                functional_t,
                OUT_TENSOR>;

        using op_type = alg::operators::linear_operator<op_impl>;
        op_type op;

        fixture() : op{
                {functional_t(SHUFFLE_TENSOR(S(1))), OUT_TENSOR(S(1))},
                {functional_t(SHUFFLE_TENSOR(TKEY(LET(1)))), OUT_TENSOR(OKEY(LET(1)))},
                {functional_t(SHUFFLE_TENSOR(TKEY(LET(2)))), OUT_TENSOR(OKEY(LET(2)))},
                {functional_t(SHUFFLE_TENSOR(TKEY(LET(3)))), OUT_TENSOR(OKEY(LET(3)))},
                {functional_t(SHUFFLE_TENSOR(TKEY(LET(4)))), OUT_TENSOR(OKEY(LET(4)))},
                {functional_t(SHUFFLE_TENSOR(TKEY(LET(5)))), []() {
                     OUT_TENSOR r;
                     for (LET l = 1; l <= 4; ++l) {
                         r += OUT_TENSOR(l, S(1));
                     }
                     return r;
                 }()}}
        {
        }
    };

    TEST_FIXTURE(fixture, test_zero)
    {
        TENSOR in;
        OUT_TENSOR expected;
        auto out = op(in);

        CHECK_EQUAL(expected, out);
    }

    TEST_FIXTURE(fixture, test_empty_word)
    {
        TENSOR in(S(1));
        OUT_TENSOR expected(S(1));

        auto out = op(in);
        CHECK_EQUAL(expected, out);
    }

    TEST_FIXTURE(fixture, test_first_letters)
    {
        for (LET l = 1; l < ALPHABET_SIZE; ++l) {
            TENSOR in(l, S(1));
            OUT_TENSOR expected(l, S(1));

            auto out = op(in);
            CHECK_EQUAL(expected, out);
        }
    }

    TEST_FIXTURE(fixture, test_final_letter)
    {
        TENSOR in(LET(ALPHABET_SIZE), S(1));
        OUT_TENSOR expected;
        for (LET l = 1; l < ALPHABET_SIZE; ++l) {
            expected += OUT_TENSOR(l, S(1));
        }
        auto out = op(in);

        CHECK_EQUAL(expected, out);
    }
}
