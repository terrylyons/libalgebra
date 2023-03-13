//
// Created by sam on 30/10/2021.
//

#include <UnitTest++.h>

#include <libalgebra/alg_types.h>
#include <libalgebra/operators.h>

#include "../../common/random_vector_generator.h"
#include "../../common/rng.h"

SUITE(tensor_left_multiplication_operator)
{

    struct fixture : alg_types<5, 5, Rational> {
        TENSOR reference_tensor;
        alg::operators::left_multiplication_operator<TENSOR> op;

        using rat_dist = la_testing::uniform_rational_distribution<S>;
        using rvg_t = la_testing::random_vector_generator<TENSOR, rat_dist>;

        std::mt19937 rng;
        rvg_t rvg;

        fixture() : reference_tensor(typename TENSOR::KEY(alg::LET(1))),
                    op(TENSOR(reference_tensor)), rng(std::random_device()()), rvg(-1, 1)
        {}
    };

    TEST_FIXTURE(fixture, test_left_multiplication_operator_action)
    {
        auto argument = rvg(rng);

        auto result = op(argument);

        CHECK_EQUAL(reference_tensor * argument, result);
    }
}

SUITE(tensor_right_multiplication_operator)
{

    struct fixture : alg_types<5, 5, Rational> {
        TENSOR reference_tensor;
        alg::operators::right_multiplication_operator<TENSOR> op;

        using rat_dist = la_testing::uniform_rational_distribution<S>;
        using rvg_t = la_testing::random_vector_generator<TENSOR, rat_dist>;

        std::mt19937 rng;
        rvg_t rvg;

        fixture() : reference_tensor(typename TENSOR::KEY(alg::LET(1))),
                    op(TENSOR(reference_tensor)), rng(std::random_device()()), rvg(-1, 1)
        {}
    };

    TEST_FIXTURE(fixture, test_right_multiplication_operator_action)
    {
        auto argument = rvg(rng);

        auto result = op(argument);

        CHECK_EQUAL(argument * reference_tensor, result);
    }
}
