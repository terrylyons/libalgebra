//
// Created by user on 04/04/23.
//
#include "width1.h"
#include <libalgebra/tensor.h>

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

    TEST_FIXTURE(Width1Tests, test_tensor_antipode)
    {
        typename TENSOR::KEY key;
        DEG d = 0;

        TENSOR tensor(key, S(d + 1));
        for (; d < DEPTH; ++d) {
            key.push_back(LET(1));
            tensor.add_scal_prod(key, S(d + 1));
        }

        auto result = antipode(tensor);

        typename TENSOR::KEY ekey;
        TENSOR expected(ekey, S(1));
        for (d = 0; d < DEPTH; ++d) {
            ekey.push_back(LET(1));
            if (d % 2 == 0) {
                expected.add_scal_prod(ekey, S(d + 1));
            }
            else {
                expected.sub_scal_prod(ekey, S(d + 1));
            }
        }
    }

    TEST_FIXTURE(Width1Tests, test_half_shuffle)
    {
        typename TENSOR::KEY key1;
        DEG d1 = 0;
        TENSOR tensor1(key1, S(d1 + 1));
        for (; d1 <DEPTH; ++d1) {
            key1.push_back(LET(1));
            tensor1.add_scal_prod(key1, S(d1 + 1));
        }

        typename TENSOR::KEY key2;
        DEG d2 = 0;
        TENSOR tensor2(key2, S(d2 + 1));
        for (; d2 < DEPTH; ++d2) {
            key2.push_back(LET(1));
            tensor2.add_scal_prod(key2, S(d2 + 1));
        }


        alg::half_shuffle_multiplication<ALPHABET_SIZE, DEPTH> mul;

        TENSOR result;
        alg::dtl::original_algebras<TENSOR, TENSOR, TENSOR> orig{ tensor1, tensor2, result};
        mul.fma(result.base_vector(),
                tensor1.base_vector(),
                tensor2.base_vector(),
                alg::mult::scalar_passthrough(),
                DEPTH,
                orig
                );

        CHECK_EQUAL(DEPTH, result.size());
    }

    TEST_FIXTURE(Width1Tests, test_area_product)
    {
        typename TENSOR::KEY key1;
        DEG d1 = 0;
        TENSOR tensor1(key1, S(d1 + 1));
        for (; d1 < DEPTH; ++d1) {
            key1.push_back(LET(1));
            tensor1.add_scal_prod(key1, S(d1 + 1));
        }

        typename TENSOR::KEY key2;
        DEG d2 = 0;
        TENSOR tensor2(key2, S(d2 + 1));
        for (; d2 < DEPTH; ++d2) {
            key2.push_back(LET(1));
            tensor2.add_scal_prod(key2, S(d2 + 1));
        }


        alg::area_tensor_multiplication <ALPHABET_SIZE, DEPTH> mul;

        TENSOR result;
        alg::dtl::original_algebras<TENSOR, TENSOR, TENSOR> orig{ tensor1, tensor2, result};
        mul.fma(result.base_vector(),
                tensor1.base_vector(),
                tensor2.base_vector(),
                alg::mult::scalar_passthrough(),
                DEPTH,
                orig
                );

        CHECK_EQUAL(0, result.size()); // ? TODO CHECK
    }
    TEST_FIXTURE(Width1Tests, test_shuffle_product)
    {
        typename SHUFFLE_TENSOR::KEY key1;
        DEG d1 = 0;
        SHUFFLE_TENSOR tensor1(key1, S(d1 + 1));
        for (; d1 < DEPTH; ++d1) {
            key1.push_back(LET(1));
            tensor1.add_scal_prod(key1, S(d1 + 1));
        }

        typename SHUFFLE_TENSOR::KEY key2;
        DEG d2 = 0;
        SHUFFLE_TENSOR tensor2(key2, S(d2 + 1));
        for (; d2 < DEPTH; ++d2) {
            key2.push_back(LET(1));
            tensor2.add_scal_prod(key2, S(d2 + 1));
        }

        auto result = tensor1*tensor2;

        CHECK_EQUAL(DEPTH+1, result.size());
    }
}
