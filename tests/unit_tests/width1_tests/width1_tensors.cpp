//
// Created by user on 04/04/23.
//
#include "width1.h"

SUITE(Width1)
{

    TEST_FIXTURE(Width1Tests, test_tensor_scalar_creation)
    {
        TENSOR t(1.0);

        typename TENSOR::KEY key;

        CHECK_EQUAL(1.0, t[key]);
    }

    TEST_FIXTURE(Width1Tests, test_tensor_letter_creation)
    {
        typename TENSOR::KEY key(LET(1));
        TENSOR t(key, 1.0);

        CHECK_EQUAL(1.0, t[key]);
    }

    TEST_FIXTURE(Width1Tests, test_tensor_addition)
    {
        TENSOR lhs(1.0);
        typename TENSOR::KEY k1(LET(1));
        TENSOR rhs(k1, 2.0);

        auto result = lhs + rhs;

        CHECK_EQUAL(1.0, result[typename TENSOR::KEY()]);
        CHECK_EQUAL(2.0, result[k1]);
    }
}
