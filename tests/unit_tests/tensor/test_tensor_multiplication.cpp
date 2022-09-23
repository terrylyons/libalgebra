//
// Created by sam on 04/02/2021.
//

#include <UnitTest++/UnitTest++.h>

#include <libalgebra/alg_types.h>
#include <libalgebra/libalgebra.h>
#include <libalgebra/polynomial_coefficients.h>

#include "../../common/time_and_details.h"

using alg::DEG;
using alg::LET;

SUITE(tensor_multiplication)
{

    struct Fixture : public alg_types<5, 5, Rational> {
        typedef alg_types<5, 5, Rational> ALG_TYPES;
        typedef typename ALG_TYPES::TENSOR TENSOR;
        typedef typename TENSOR::BASIS TBASIS;
        typedef typename TENSOR::KEY KEY;

        const TENSOR tunit;
        const TENSOR tzero;

        Fixture() : tunit(KEY()), tzero()
        {}

        KEY make_key(const LET* arg, const std::size_t N)
        {
            KEY k;
            for (std::size_t i = 0; i < N; ++i) {
                k.push_back(arg[i]);
            }
            return k;
        }
    };

#define ADD_KEY(N, ...)                       \
    {                                         \
        LET tmp[N] = {__VA_ARGS__};           \
        expected += TENSOR(make_key(tmp, N)); \
    }

    TEST_FIXTURE(Fixture, test_product_tunit_zero)
    {
        TEST_DETAILS();
        TENSOR lhs(tunit), rhs, expected;

        CHECK_EQUAL(expected, lhs * rhs);
    }

    TEST_FIXTURE(Fixture, test_product_zero_tunit)
    {
        TEST_DETAILS();
        TENSOR lhs, rhs(tunit), expected;

        CHECK_EQUAL(expected, lhs * rhs);
    }

    TEST_FIXTURE(Fixture, test_product_tunit_tunit)
    {
        TEST_DETAILS();

        TENSOR lhs(tunit), rhs(tunit);

        CHECK_EQUAL(tunit, lhs * rhs);
    }

    TEST_FIXTURE(Fixture, test_product_unidim_deg_1_tunit)
    {
        TEST_DETAILS();

        LET k1[] = {1};

        TENSOR lhs(tunit), rhs(make_key(k1, 1));

        TENSOR expected(rhs);
        CHECK_EQUAL(expected, lhs * rhs);
    }

    TEST_FIXTURE(Fixture, test_product_unidim_deg_1_1)
    {
        TEST_DETAILS();

        LET k1[] = {1};
        LET k2[] = {2};
        TENSOR lhs(make_key(k1, 1)), rhs(make_key(k2, 1));

        LET ek[] = {1, 2};
        TENSOR expected(make_key(ek, 2));

        CHECK_EQUAL(expected, lhs * rhs);
    }

    TEST_FIXTURE(Fixture, test_product_unidim_deg_1_deg_2)
    {
        TEST_DETAILS();

        LET k1[] = {1};
        LET k2[] = {2, 3};
        TENSOR lhs(make_key(k1, 1)), rhs(make_key(k2, 2));

        LET ek[] = {1, 2, 3};
        TENSOR expected(make_key(ek, 3));

        CHECK_EQUAL(expected, lhs * rhs);
    }

    TEST_FIXTURE(Fixture, test_product_unidim_deg_2_deg_1)
    {
        TEST_DETAILS();

        LET k1[] = {1, 2};
        LET k2[] = {3};
        TENSOR lhs(make_key(k1, 2)), rhs(make_key(k2, 1));

        LET ek[] = {1, 2, 3};
        TENSOR expected(make_key(ek, 3));

        CHECK_EQUAL(expected, lhs * rhs);
    }

    TEST_FIXTURE(Fixture, test_product_unidim_deg_3_deg_3_overflow_removed)
    {
        TEST_DETAILS();

        LET k1[] = {1, 2, 3};
        LET k2[] = {1, 2, 3};
        TENSOR lhs(make_key(k1, 3)), rhs(make_key(k2, 3));

        TENSOR expected(tzero);

        CHECK_EQUAL(expected, lhs * rhs);
    }

    TEST_FIXTURE(Fixture, test_product_multiple_terms)
    {
        TEST_DETAILS();

        TENSOR lhs(tunit), rhs(tunit);

        LET letters[] = {1, 2, 3, 4, 4};
        // 1() + 1(1) + 1(1,2) + 1(1,2,3) + 1(1,2,3,4)
        for (DEG d = 1; d < 5; ++d) {
            lhs += TENSOR(make_key(&letters[0], d));
        }

        // 1() + 1(4) + 1(3,4) + 1(2,3,4) + 1(1,2,3,4)
        for (DEG d = 1; d < 5; ++d) {
            rhs += TENSOR(make_key(&letters[4 - d], d));
        }

        TENSOR expected(tunit);

        const S one(1);

        // DEG 5
        ADD_KEY(5, 1, 2, 3, 4, 4);
        ADD_KEY(5, 1, 2, 3, 3, 4);
        ADD_KEY(5, 1, 2, 2, 3, 4);
        ADD_KEY(5, 1, 1, 2, 3, 4);

        /*
    expected += TENSOR(make_key((LET[]) {1,2,3,4,4}, 5));
    expected += TENSOR(make_key((LET[]) {1,2,3,3,4}, 5));
    expected += TENSOR(make_key((LET[]) {1,2,2,3,4}, 5));
    expected += TENSOR(make_key((LET[]) {1,1,2,3,4}, 5));
    */
        // DEG 4
        // (1,2,3,4) = () * (1,2,3,4), (1) * (2,3,4), (1,2) * (3,4), (1,2,3) * (4), (1,2,3,4) * ()
        {
            LET tmp[] = {1, 2, 3, 4};
            expected.add_scal_prod(make_key(tmp, 4), S(5));
        }

        // DEG 3
        ADD_KEY(3, 1, 2, 3);
        ADD_KEY(3, 1, 2, 4);
        ADD_KEY(3, 1, 3, 4);
        ADD_KEY(3, 2, 3, 4);
        /*
    expected += TENSOR(make_key((LET[]) {1,2,3}, 3));
    expected += TENSOR(make_key((LET[]) {1,2,4}, 3));
    expected += TENSOR(make_key((LET[]) {1,3,4}, 3));
    expected += TENSOR(make_key((LET[]) {2,3,4}, 3));
*/

        // DEG 2
        ADD_KEY(2, 1, 2);
        ADD_KEY(2, 1, 4);
        ADD_KEY(2, 3, 4);
        /*
    expected += TENSOR(make_key((LET[]) {1,2}, 2));
    expected += TENSOR(make_key((LET[]) {1,4}, 2));
    expected += TENSOR(make_key((LET[]) {3,4}, 2));
*/

        // DEG 1
        expected += TENSOR(KEY(LET(1)));
        expected += TENSOR(KEY(LET(4)));

        CHECK_EQUAL(expected, lhs * rhs);
    }

#undef ADD_KEY

    struct PolyMultiplicationTests {
        static constexpr DEG width = 5;
        static constexpr DEG depth = 5;

        using coeff_ring = alg::coefficients::rational_poly_ring;
        using rational_type = typename coeff_ring::Q;
        using scalar_type = typename coeff_ring::S;

        using poly_basis = alg::poly_basis;
        using poly_key = typename poly_basis::KEY;

        using tensor_basis = alg::tensor_basis<width, depth>;
        using traditional_multiplication = alg::traditional_free_tensor_multiplication<width, depth>;
        using tiled_multiplication = alg::tiled_free_tensor_multiplication<width, depth>;
        using tiled2_multiplication = alg::tiled_free_tensor_multiplication<width, depth, 2>;

        using key_type = typename tensor_basis::KEY;
        tensor_basis basis;

        using sparse_traditional_tensor = alg::algebra<tensor_basis, coeff_ring, traditional_multiplication, alg::vectors::sparse_vector>;
        using dense_traditional_tensor = alg::algebra<tensor_basis, coeff_ring, traditional_multiplication, alg::vectors::dense_vector>;

        scalar_type construct_expected(key_type key, LET lhs_offset, LET rhs_offset) const
        {
            auto deg = key.size();

            scalar_type result;
            for (DEG i = 0; i <= deg; ++i) {
                auto right(key);
                auto left = right.split_n(i);
                result += to_poly_key(left, lhs_offset) * to_poly_key(right, rhs_offset);
            }
            return result;
        }

        scalar_type to_poly_key(key_type key, LET offset) const
        {
            LET count = 0;
            while (key.size() > 0) {
                count *= 10;
                count += key.FirstLetter();
                key = key.rparent();
            }
            return scalar_type(count + offset, 1);
        }

        template<typename Mul, template<typename, typename, typename...> class VT>
        alg::algebra<tensor_basis, coeff_ring, Mul, VT> generic_tensor(LET offset = 0) const
        {
            alg::algebra<tensor_basis, coeff_ring, Mul, VT> result;

            for (auto key : basis.iterate_keys()) {
                result[key] = to_poly_key(key, offset);
            }
            return result;
        }
    };

#define LA_TESTING_TENSOR_MUL_MAKE_TESTS(MUL, VEC)                                            \
    TEST_FIXTURE(PolyMultiplicationTests, VEC##_##MUL##_fma)                                  \
    {                                                                                         \
        auto lhs = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>(1000000); \
        auto rhs = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>(2000000); \
                                                                                              \
        auto result = lhs * rhs;                                                              \
        REQUIRE CHECK_EQUAL(3906, result.size());                                             \
        for (auto item : result) {                                                            \
            CHECK_EQUAL(construct_expected(item.key(), 1000000, 2000000), item.value());      \
        }                                                                                     \
    }                                                                                         \
                                                                                              \
    TEST_FIXTURE(PolyMultiplicationTests, VEC##_##MUL##_inplace_mul)                          \
    {                                                                                         \
        auto lhs = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>(1000000); \
        auto rhs = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>(2000000); \
                                                                                              \
        lhs *= rhs;                                                                           \
        REQUIRE CHECK_EQUAL(3906, lhs.size());                                                \
        for (auto item : lhs) {                                                               \
            CHECK_EQUAL(construct_expected(item.key(), 1000000, 2000000), item.value());      \
        }                                                                                     \
    }

    LA_TESTING_TENSOR_MUL_MAKE_TESTS(traditional, dense)
    LA_TESTING_TENSOR_MUL_MAKE_TESTS(tiled, dense)
    LA_TESTING_TENSOR_MUL_MAKE_TESTS(tiled2, dense)
}
