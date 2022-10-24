//
// Created by sam on 24/01/2022.
//

#ifndef LIBALGEBRA_UNITAL_COEFFICIENT_TEST_SUITE_H
#define LIBALGEBRA_UNITAL_COEFFICIENT_TEST_SUITE_H

#include <UnitTest++.h>
#include <libalgebra/coefficients.h>
#include <libalgebra/libalgebra.h>
#include <multi_test.h>

using alg::DEG;

NEW_AUTO_SUITE(unital_coefficients_tests, typename Coeffs)
{

    ADD_TEST(test_one)
    {
        typename Coeffs::S expected(1);
        CHECK_EQUAL(expected, Coeffs::one);
    };

    ADD_TEST(test_mone)
    {
        typename Coeffs::S expected(-1);
        CHECK_EQUAL(expected, Coeffs::mone);
    };

    ADD_TEST(test_one_plus_mone)
    {
        auto expected = Coeffs::zero;
        CHECK_EQUAL(expected, Coeffs::add(Coeffs::one, Coeffs::mone));
    };

    ADD_TEST(test_one_minus_one)
    {
        CHECK_EQUAL(Coeffs::mone, Coeffs::uminus(Coeffs::one));
    };

    ADD_TEST(test_one_subtract_one)
    {
        CHECK_EQUAL(Coeffs::zero, Coeffs::sub(Coeffs::one, Coeffs::one));
    };

    ADD_TEST(test_one_mul_zero)
    {
        CHECK_EQUAL(Coeffs::zero, Coeffs::mul(Coeffs::one, Coeffs::zero));
    };

    ADD_TEST(test_one_mul_one)
    {
        CHECK_EQUAL(Coeffs::one, Coeffs::mul(Coeffs::one, Coeffs::one));
    };

    ADD_TEST(test_one_mul_mone)
    {
        CHECK_EQUAL(Coeffs::mone, Coeffs::mul(Coeffs::one, Coeffs::mone));
    };

    ADD_TEST(test_one_div_one)
    {
        CHECK_EQUAL(Coeffs::one, Coeffs::div(Coeffs::one, Coeffs::one));
    };

    ADD_TEST(test_one_div_mone)
    {
        CHECK_EQUAL(Coeffs::mone, Coeffs::div(Coeffs::one, Coeffs::mone));
    };

    ADD_TEST(test_constructible_from_degree_type)
    {
        typename Coeffs::S expected(2);

        CHECK_EQUAL(expected, Coeffs::from(DEG(2)));
    };
}

#endif//LIBALGEBRA_UNITAL_COEFFICIENT_TEST_SUITE_H
