//
// Created by sam on 24/01/2022.
//

#ifndef LIBALGEBRA_STANDARD_COEFFICIENT_TEST_SUITE_H
#define LIBALGEBRA_STANDARD_COEFFICIENT_TEST_SUITE_H

#include <UnitTest++.h>
#include <multi_test.h>
#include <random_coeffs.h>

NEW_AUTO_SUITE(standard_coefficient_tests, typename Coeffs)
{

    ADD_TEST(test_zero)
    {
        typename Coeffs::S expected(0);

        CHECK_EQUAL(expected, Coeffs::zero);
    };

    ADD_TEST(test_add_uminus) {
        la_testing::random_coeff_generator<Coeffs> dist(-1, 1);
        std::mt19937 rng(std::random_device{}());

        auto x = dist(rng), y = dist(rng);

        CHECK_EQUAL(
                Coeffs::sub(x, y),
                Coeffs::add(x, Coeffs::uminus(y))
                );
    };


}

#endif//LIBALGEBRA_STANDARD_COEFFICIENT_TEST_SUITE_H
