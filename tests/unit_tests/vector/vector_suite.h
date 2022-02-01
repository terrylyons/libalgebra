//
// Created by sam on 01/02/2022.
//

#ifndef LIBALGEBRA_VECTOR_SUITE_H
#define LIBALGEBRA_VECTOR_SUITE_H

#include <libalgebra/coefficients/coefficients.h>
#include <libalgebra/coefficients/rational_coefficients.h>
#include <libalgebra/implementation_types.h>
#include <libalgebra/tensor_basis.h>

#include <UnitTest++/UnitTest++.h>

#include <random>

#include <multi_test.h>
#include <random_vector_generator.h>
#include <rng.h>

using alg::LET;

NEW_AUTO_SUITE(vector_suite, template<typename, typename, typename...> class VectorType)
{
    using basis_type = alg::tensor_basis<5, 2>;
    using coeff_type = alg::coefficients::rational_field;
    using tensor_vec = VectorType<basis_type, coeff_type>;

    using tkey_type = typename basis_type::KEY;

    using S = typename alg::coefficients::rational_field::SCA;

    using skip_dist = std::uniform_int_distribution<alg::DIMN>;
    using random_vect_gen = la_testing::random_vector_generator<tensor_vec,
                                                                la_testing::uniform_rational_distribution<S>, skip_dist>;

    ADD_TEST(test_default_constructor_equal)
    {
        tensor_vec left, right;

        CHECK_EQUAL(left, right);
    };

    ADD_TEST(test_size_empty_vector)
    {
        tensor_vec arg;

        CHECK_EQUAL(0, arg.size());
    };

    ADD_TEST(test_size_one_emptyword)
    {
        tensor_vec arg(tkey_type(), S(1));

        CHECK_EQUAL(1, arg.size());
    };

    ADD_TEST(test_size_two_letters)
    {
        tensor_vec arg;
        arg.add_scal_prod(tkey_type(LET(1)), S(1));
        arg.add_scal_prod(tkey_type(LET(2)), S(2));

        CHECK_EQUAL(2, arg.size());
    };

    ADD_TEST(test_degree_empty_vector)
    {
        tensor_vec arg;

        CHECK_EQUAL(0, arg.degree());
    };

    ADD_TEST(test_degree_emptyword)
    {
        tensor_vec arg(tkey_type(), S(1));

        CHECK_EQUAL(0, arg.degree());
    };

    ADD_TEST(test_degree_letter)
    {
        tensor_vec arg(tkey_type(LET(1)), S(1));

        CHECK_EQUAL(1, arg.degree());
    };

    ADD_TEST(test_degree_length_2_word)
    {
        tkey_type key{LET(1), LET(2)};
        tensor_vec arg(key, S(1));

        CHECK_EQUAL(2, arg.degree());
    };

    ADD_TEST(test_degree_multiple_words)
    {
        tensor_vec arg(tkey_type(), S(1));
        arg.add_scal_prod(tkey_type(LET(1)), S(2));
        arg.add_scal_prod(tkey_type{LET(1), LET(2)}, S(3));

        CHECK_EQUAL(2, arg.degree());
    };

    ADD_TEST(test_addition_neutral_neutral)
    {
        tensor_vec left, right, expected;

        CHECK_EQUAL(expected, left + right);
    };

    ADD_TEST(test_addition_neutral_random)
    {
        std::mt19937 rng(std::random_device{}());
        random_vect_gen gen(10UL, -2, 2);
        tensor_vec left, right(gen(rng)), expected(right);

        auto result = left + right;
        CHECK_EQUAL(expected, result);
    };

    ADD_TEST(test_inplace_vs_external_addition)
    {
        std::mt19937 rng(std::random_device{}());
        random_vect_gen gen(10UL, -2, 2);

        tensor_vec left_inplace(gen(rng)), left(left_inplace), right(gen(rng));

        left_inplace += right;
        CHECK_EQUAL(left_inplace, left + right);
    };

    ADD_TEST(test_subtraction_neutral_neutral)
    {
        tensor_vec left, right, expected;

        CHECK_EQUAL(expected, left - right);
    };

    ADD_TEST(test_subtraction_neutral_random_vs_uminus)
    {
        std::mt19937 rng(std::random_device{}());
        random_vect_gen gen(10UL, -2, 2);
        tensor_vec left, right(gen(rng)), expected(-right);

        CHECK_EQUAL(expected, left - right);
    };

    ADD_TEST(test_inplace_vs_external_subtraction)
    {
        std::mt19937 rng(std::random_device{}());
        random_vect_gen gen(10UL, -2, 2);

        tensor_vec left_inplace(gen(rng)), left(left_inplace), right(gen(rng));

        left_inplace -= right;
        CHECK_EQUAL(left_inplace, left - right);
    };

    ADD_TEST(test_fused_add_scal_prod_vs_separate)
    {
        std::mt19937 rng(std::random_device{}());
        random_vect_gen gen(10UL, -2, 2);

        tensor_vec left1(gen(rng)), left2(left1), right(gen(rng));
        S scalar(2);

        left1.add_scal_prod(right, scalar);

        left2 += (right * scalar);

        CHECK_EQUAL(left1, left2);
    };

    ADD_TEST(test_fused_sub_scal_prod_vs_separate)
    {
        std::mt19937 rng(std::random_device{}());
        random_vect_gen gen(10UL, -2, 2);

        tensor_vec left1(gen(rng)), left2(left1), right(gen(rng));
        S scalar(2);

        left1.sub_scal_prod(right, scalar);

        left2 -= (right * scalar);

        CHECK_EQUAL(left1, left2);
    };

    ADD_TEST(test_fused_add_scal_div_vs_separate)
    {
        std::mt19937 rng(std::random_device{}());
        random_vect_gen gen(10UL, -2, 2);

        tensor_vec left1(gen(rng)), left2(left1), right(gen(rng));
        S scalar(2);

        left1.add_scal_div(right, scalar);

        left2 += (right / scalar);

        CHECK_EQUAL(left1, left2);
    };

    ADD_TEST(test_fused_sub_scal_div_vs_separate)
    {
        std::mt19937 rng(std::random_device{}());
        random_vect_gen gen(10UL, -2, 2);

        tensor_vec left1(gen(rng)), left2(left1), right(gen(rng));
        S scalar(2);

        left1.sub_scal_div(right, scalar);

        left2 -= (right / scalar);

        CHECK_EQUAL(left1, left2);
    };
}

#endif//LIBALGEBRA_VECTOR_SUITE_H
