//
// Created by sam on 30/10/2021.
//

#include <UnitTest++/UnitTest++.h>

#include <libalgebra/alg_types.h>
#include <libalgebra/operators.h>

#include "../../common/random_vector_generator.h"
#include "../../common/rng.h"

SUITE(tensor_left_multiplication_operator) {

    struct fixture : alg_types<5, 5, Rational>
    {
        TENSOR reference_tensor;
        alg::operators::left_multiplication_operator<TENSOR> op;

        using rat_dist = la_unittests::
        using rvg = la_unittests::random_vector_generator<TENSOR>

        fixture() : reference_tensor(typename TENSOR::KEY(alg::LET(1))), op(TENSOR(reference_tensor))
        {}
    };

    TEST_FIXTURE(fixture, test_left_multiplication_operator_action)
    {
        auto argument =
    }

}