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

using alg::LET;

NEW_AUTO_SUITE(vector_suite, template<typename, typename, typename...> class VectorType)
{
    using basis_type = alg::tensor_basis<5, 5>;
    using coeff_type = alg::coefficients::rational_field;
    using tensor_vec = VectorType<basis_type, coeff_type>;

    using tkey_type = typename basis_type::KEY;

    std::mt19937 rng(std::random_device{}());

    using S = typename alg::coefficients::rational_field::SCA;

    la_testing::random_vector_generator<tensor_vec> vector_generator;

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
}

#endif//LIBALGEBRA_VECTOR_SUITE_H
