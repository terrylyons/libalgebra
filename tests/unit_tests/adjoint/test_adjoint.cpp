#include <UnitTest++.h>

#include <iostream>
#include <vector>

#include <libalgebra/libalgebra.h>
#include <libalgebra/tensor.h>
#include "../../common/random_vector_generator.h"
#include "../../common/rng.h"
#include <libalgebra/polynomial_coefficients.h>
#include <libalgebra/operators.h>

using namespace alg;

SUITE(Adjoint)
{
    template<DEG Width, DEG Depth>
    struct PolyMultiplicationTests {
        static constexpr DEG width = Width;
        static constexpr DEG depth = Depth;

        using coeff_ring = alg::coefficients::rational_poly_ring;
        using rational_type = typename coeff_ring::Q;
        using scalar_type = typename coeff_ring::S;

        using poly_basis = alg::poly_basis;
        using poly_key = typename poly_basis::KEY;

        using tensor_basis = alg::free_tensor_basis<width, depth>;
        using traditional_multiplication = alg::traditional_free_tensor_multiplication<width, depth>;
        using tiled_multiplication = alg::tiled_free_tensor_multiplication<width, depth>;
        using tiled2_multiplication = alg::tiled_free_tensor_multiplication<width, depth, 2>;

        using key_type = typename tensor_basis::KEY;
        tensor_basis basis;

        using tensor_type = alg::free_tensor<coeff_ring, width, depth>;

        using shuffle_tensor_type = alg::shuffle_tensor<coeff_ring, width, depth>;

        template <typename T>
        using operator_type = operators::dtl::adjoint_of_left_multiplication_operator_impl<T>;

        using sparse_traditional_tensor = alg::algebra<tensor_basis, coeff_ring, traditional_multiplication, alg::vectors::sparse_vector>;
        using dense_traditional_tensor = alg::algebra<tensor_basis, coeff_ring, traditional_multiplication, alg::vectors::dense_vector>;

        scalar_type construct_expected(key_type key, LET lhs_offset, LET rhs_offset, DEG lhs_max = depth, DEG rhs_max = depth) const
        {
            auto deg = key.size();
            auto dmin = (deg >= rhs_max) ? deg - rhs_max : DEG(0);
            auto dmax = std::min(deg, lhs_max);

            scalar_type result;
            for (DEG i = dmin; i <= dmax; ++i) {
                auto right(key);
                auto left = right.split_n(i);
                result += to_poly_key(left, lhs_offset) * to_poly_key(right, rhs_offset);
            }
            return result;
        }

        scalar_type construct_expected_rzu(key_type key, LET lhs_offset, LET rhs_offset) const
        {
            auto deg = key.size();

            scalar_type result;
            for (DEG i = 0; i < deg; ++i) {
                auto right(key);
                auto left = right.split_n(i);
                result += to_poly_key(left, lhs_offset) * to_poly_key(right, rhs_offset);
            }
            return result;
        }
        scalar_type construct_expected_lzu(key_type key, LET lhs_offset, LET rhs_offset) const
        {
            auto deg = key.size();

            scalar_type result;
            for (DEG i = 1; i <= deg; ++i) {
                auto right(key);
                auto left = right.split_n(i);
                result += to_poly_key(left, lhs_offset) * to_poly_key(right, rhs_offset);
            }
            return result;
        }

        scalar_type to_poly_key(key_type key, LET offset) const
        {
            constexpr LET digits_per_letter = alg::integer_maths::logN(LET(width), LET(10)) + 1;
            constexpr LET letter_offset = alg::power(LET(10), digits_per_letter);

            auto offset_digits = alg::integer_maths::logN(LET(offset), LET(10)) + LET(1);
            if (offset_digits <= depth * digits_per_letter) {
                offset *= alg::power(LET(10), depth * digits_per_letter - offset_digits + 1);
            }

            LET count = 0;
            while (key.size() > 0) {
                count *= letter_offset;
                count += key.FirstLetter();
                key = key.rparent();
            }
            return scalar_type(count + offset, 1);
        }

        template<typename Mul, template<typename, typename, typename...> class VT>
        alg::algebra<tensor_basis, coeff_ring, Mul, VT> generic_tensor(LET offset = 0, DEG max_degree = depth) const
        {
            alg::algebra<tensor_basis, coeff_ring, Mul, VT> result;

            for (auto key : basis.iterate_keys_to_deg(max_degree + 1)) {
                result[key] = to_poly_key(key, offset);
            }
            return result;
        }

        shuffle_tensor_type generic_shuffle_tensor(LET offset = 0, DEG max_degree = depth) const
        {
            shuffle_tensor_type result;

            for (auto key : basis.iterate_keys_to_deg(max_degree + 1)) {
                result[key] = to_poly_key(key, offset);
            }
            return result;
        }

        template <template <typename, typename, typename...> class VT>
        alg::shuffle_tensor<coeff_ring, Width, Depth, VT> generic_shuffle_tensor(LET offset = 0, DEG max_degree = depth) const
        {
            alg::shuffle_tensor<coeff_ring, Width, Depth, VT> result;

            for (auto key : basis.iterate_keys_to_deg(max_degree + 1)) {
                result[key] = to_poly_key(key, offset);
            }
            return result;
        }

        tensor_type generic_tensor(LET offset = 0, DEG max_degree = depth) const
        {
            tensor_type result;

            for (auto key : basis.iterate_keys_to_deg(max_degree + 1)) {
                result[key] = to_poly_key(key, offset);
            }
            return result;
        }

        alg::free_tensor<coeff_ring, width, depth, alg::vectors::dense_vector>
        generic_d_free_tensor(LET offset = 0, DEG max_degree = depth)
        {
            alg::free_tensor<coeff_ring, width, depth, alg::vectors::dense_vector> result;

            for (auto key : basis.iterate_keys_to_deg(max_degree + 1)) {
                result[key] = to_poly_key(key, offset);
            }
            return result;
        }

        template <typename Shuffle=shuffle_tensor_type>
        Shuffle construct_expected(LET tensor_offset, LET shuffle_offset) {
            Shuffle result;

            for (auto&& key : basis.iterate_keys()) {
                for (DEG i=0; i<=key.size(); ++i) {
                    auto right(key);
                    auto left = right.split_n(i);
                    result.add_scal_prod(right, to_poly_key(left, tensor_offset) * to_poly_key(key, shuffle_offset));
                }
            }

            return result;
        }
        shuffle_tensor_type construct_expected(LET op_offset, LET arg_offset, std::initializer_list<key_type> keys) const
        {
            shuffle_tensor_type result;

            for (auto&& key : basis.iterate_keys()) {
                for (const auto& op_key : keys) {
                    auto right(key);
                    auto left = right.split_n(op_key.size());
                    if (left == op_key) {
                        result.add_scal_prod(right, to_poly_key(op_key, op_offset) * to_poly_key(arg_offset, key));
                    }
                }
            }
            return result;
        }
    };

    using PolyMultiplicationTests55 = PolyMultiplicationTests<5, 5>;
    TEST_FIXTURE(PolyMultiplicationTests55, IdentityTest)
    {

        tensor_type this_tensor(scalar_type(1));
        shuffle_tensor_type this_shuffle_tensor(scalar_type(1));

        operator_type<tensor_type> adjoint(this_tensor);

        CHECK_EQUAL(this_shuffle_tensor, adjoint(this_shuffle_tensor));


    }

    TEST_FIXTURE(PolyMultiplicationTests55, PolynomialTest)
    {

        tensor_type this_tensor(scalar_type(1));
        shuffle_tensor_type this_shuffle_tensor = generic_shuffle_tensor(2);

        operator_type<tensor_type> adjoint(this_tensor);

        CHECK_EQUAL(this_shuffle_tensor, adjoint(this_shuffle_tensor));


    }


    TEST_FIXTURE(PolyMultiplicationTests55, FullTest) {
        tensor_type this_tensor(generic_tensor(1));
        shuffle_tensor_type  this_shuffle_tensor(generic_shuffle_tensor(2));

        operator_type<tensor_type> adjoint(this_tensor);

        auto expected = construct_expected(1, 2);

        CHECK_EQUAL(expected, adjoint(this_shuffle_tensor));
    }

    TEST_FIXTURE(PolyMultiplicationTests55, FullTestDense) {
        auto this_tensor(generic_tensor<traditional_multiplication, alg::vectors::dense_vector>(1));
        auto this_shuffle_tensor(generic_shuffle_tensor<alg::vectors::dense_vector>(2));

        operator_type<decltype(this_tensor)> adjoint(this_tensor);

        auto expected = construct_expected<decltype(this_shuffle_tensor)>(1, 2);

        auto adj = adjoint(this_shuffle_tensor);
//        CHECK_EQUAL(expected, adj);
        CHECK_EQUAL(decltype(this_shuffle_tensor)(), expected-adj);
    }

}// SUITE adjoint
