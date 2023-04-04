//
// Created by sam on 04/02/2021.
//

#include <UnitTest++.h>

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

        template<typename M, template<typename, typename, typename...> class V>
        using tensor_type = alg::algebra<tensor_basis, coeff_ring, M, V>;

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

        alg::free_tensor<coeff_ring, width, depth, alg::vectors::dense_vector>
        generic_d_free_tensor(LET offset = 0, DEG max_degree = depth)
        {
            alg::free_tensor<coeff_ring, width, depth, alg::vectors::dense_vector> result;

            for (auto key : basis.iterate_keys_to_deg(max_degree + 1)) {
                result[key] = to_poly_key(key, offset);
            }
            return result;
        }
    };

    using PolyMultiplicationTests55 = PolyMultiplicationTests<5, 5>;

#define LA_TESTING_TENSOR_MAX_DEG_FMA(MUL, VEC, DEGREE)                                              \
    TEST_FIXTURE(PolyMultiplicationTests55, VEC##_##MUL##_fma_max_deg_##DEGREE)                      \
    {                                                                                                \
        using mtraits = alg::dtl::multiplication_traits<MUL##_multiplication>;                       \
        auto lhs = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>(1000000);        \
        auto rhs = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>(2000000);        \
                                                                                                     \
        const auto& mul = lhs.multiplication();                                                      \
        tensor_type<MUL##_multiplication, alg::vectors::VEC##_vector> result;                        \
                                                                                                     \
        mtraits::multiply_and_add(mul, result, lhs, rhs, alg::mult::scalar_passthrough(), DEGREE);   \
                                                                                                     \
        CHECK_EQUAL(basis.start_of_degree(DEGREE + 1), result.size());                               \
        for (auto item : result) {                                                                   \
            if (basis.degree(item.key()) < DEGREE + 1) {                                             \
                REQUIRE CHECK_EQUAL(construct_expected(item.key(), 1000000, 2000000), item.value()); \
            }                                                                                        \
            else {                                                                                   \
                REQUIRE CHECK_EQUAL(scalar_type(), item.value());                                    \
            }                                                                                        \
        }                                                                                            \
    }

#define LA_TESTING_TENSOR_FMA_LOWER_DEG(MUL, VEC, DEGREE, LDEG, RDEG)                                        \
    TEST_FIXTURE(PolyMultiplicationTests55, VEC##_##MUL##_fma_lhs_smaller_##LDEG##_##RDEG)                   \
    {                                                                                                        \
        using mtraits = alg::dtl::multiplication_traits<MUL##_multiplication>;                               \
        auto lhs = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>(1000000, LDEG);          \
        auto rhs = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>(2000000, RDEG);          \
        REQUIRE CHECK_EQUAL(tensor_basis::start_of_degree(LDEG + 1), lhs.size());                            \
        REQUIRE CHECK_EQUAL(tensor_basis::start_of_degree(RDEG + 1), rhs.size());                            \
        const auto& mul = lhs.multiplication();                                                              \
        tensor_type<MUL##_multiplication, alg::vectors::VEC##_vector> result;                                \
                                                                                                             \
        mtraits::multiply_and_add(mul, result, lhs, rhs, alg::mult::scalar_passthrough());                   \
                                                                                                             \
        CHECK_EQUAL(basis.start_of_degree(std::min(LDEG + RDEG, DEGREE) + 1), result.size());                \
        for (auto item : result) {                                                                           \
            REQUIRE CHECK_EQUAL(construct_expected(item.key(), 1000000, 2000000, LDEG, RDEG), item.value()); \
        }                                                                                                    \
    }

#define LA_TESTING_TENSOR_MAX_DEG_MUL_INPLACE(MUL, VEC, DEGREE)                                      \
    TEST_FIXTURE(PolyMultiplicationTests55, VEC##_##MUL##_minplace_max_deg_##DEGREE)                 \
    {                                                                                                \
        using mtraits = alg::dtl::multiplication_traits<MUL##_multiplication>;                       \
        auto lhs = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>(1000000);        \
        auto rhs = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>(2000000);        \
                                                                                                     \
        const auto& mul = lhs.multiplication();                                                      \
                                                                                                     \
        mtraits::multiply_inplace(mul, lhs, rhs, alg::mult::scalar_passthrough(), DEGREE);           \
                                                                                                     \
        CHECK_EQUAL(3906, lhs.size());                                                               \
        for (auto item : lhs) {                                                                      \
            if (basis.degree(item.key()) < (DEGREE) + 1) {                                           \
                REQUIRE CHECK_EQUAL(construct_expected(item.key(), 1000000, 2000000), item.value()); \
            }                                                                                        \
            else {                                                                                   \
                REQUIRE CHECK_EQUAL(to_poly_key(item.key(), 1000000), item.value());                 \
            }                                                                                        \
        }                                                                                            \
    }

#define LA_TESTING_TENSOR_MAX_DEG_MSD(MUL, VEC, DEGREE)                                                      \
    TEST_FIXTURE(PolyMultiplicationTests55, VEC##_##MUL##_msd_max_deg_rzu_##DEGREE)                          \
    {                                                                                                        \
        auto lhs = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>(1000000);                \
        auto rhs = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>(2000000);                \
        rhs[key_type()] = scalar_type();                                                                     \
                                                                                                             \
        lhs.mul_scal_div(rhs, rational_type(2), DEGREE);                                                     \
                                                                                                             \
        CHECK_EQUAL(3905, lhs.size());                                                                       \
        for (auto item : lhs) {                                                                              \
            if (item.key() == key_type()) {                                                                  \
                REQUIRE CHECK_EQUAL(scalar_type(), item.value());                                            \
            }                                                                                                \
            else if (basis.degree(item.key()) < (DEGREE) + 1) {                                              \
                REQUIRE CHECK_EQUAL(construct_expected_rzu(item.key(), 1000000, 2000000) / 2, item.value()); \
            }                                                                                                \
            else {                                                                                           \
                REQUIRE CHECK_EQUAL(to_poly_key(item.key(), 1000000), item.value());                         \
            }                                                                                                \
        }                                                                                                    \
    }

#define LA_TESTING_TENSOR_MUL_MAKE_TESTS(MUL, VEC)                                                        \
    TEST_FIXTURE(PolyMultiplicationTests55, VEC##_##MUL##_fma)                                            \
    {                                                                                                     \
        auto lhs = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>(1000000);             \
        auto rhs = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>(2000000);             \
                                                                                                          \
        auto result = lhs * rhs;                                                                          \
        CHECK_EQUAL(3906, result.size());                                                         \
        for (auto item : result) {                                                                        \
            REQUIRE CHECK_EQUAL(construct_expected(item.key(), 1000000, 2000000), item.value());          \
        }                                                                                                 \
    }                                                                                                     \
                                                                                                          \
    TEST_FIXTURE(PolyMultiplicationTests55, VEC##_##MUL##_inplace_mul)                                    \
    {                                                                                                     \
        auto lhs = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>(1000000);             \
        auto rhs = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>(2000000);             \
                                                                                                          \
        lhs *= rhs;                                                                                       \
        CHECK_EQUAL(3906, lhs.size());                                                            \
        for (auto item : lhs) {                                                                           \
            REQUIRE CHECK_EQUAL(construct_expected(item.key(), 1000000, 2000000), item.value());          \
        }                                                                                                 \
    }                                                                                                     \
                                                                                                          \
    TEST_FIXTURE(PolyMultiplicationTests55, VEC##_##MUL##_inplace_mul_rhs_zero_unit)                      \
    {                                                                                                     \
        auto lhs = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>(1000000);             \
        auto rhs = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>(2000000);             \
        rhs[key_type()] = scalar_type();                                                                  \
                                                                                                          \
        lhs *= rhs;                                                                                       \
        CHECK_EQUAL(3905, lhs.size());                                                            \
        for (auto item : lhs) {                                                                           \
            if (basis.degree(item.key()) == 0) {                                                          \
                REQUIRE CHECK_EQUAL(scalar_type(), item.value());                                         \
            }                                                                                             \
            else {                                                                                        \
                REQUIRE CHECK_EQUAL(construct_expected_rzu(item.key(), 1000000, 2000000), item.value());  \
            }                                                                                             \
        }                                                                                                 \
    }                                                                                                     \
                                                                                                          \
    TEST_FIXTURE(PolyMultiplicationTests55, VEC##_##MUL##_inplace_mul_lhs_zero_unit)                      \
    {                                                                                                     \
        auto lhs = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>(1000000);             \
        auto rhs = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>(2000000);             \
        lhs[key_type()] = scalar_type();                                                                  \
                                                                                                          \
        lhs *= rhs;                                                                                       \
        CHECK_EQUAL(3905, lhs.size());                                                            \
        for (auto item : lhs) {                                                                           \
            if (basis.degree(item.key()) == 0) {                                                          \
                REQUIRE CHECK_EQUAL(scalar_type(), item.value());                                                 \
            }                                                                                             \
            else {                                                                                        \
                REQUIRE CHECK_EQUAL(construct_expected_lzu(item.key(), 1000000, 2000000), item.value());          \
            }                                                                                             \
        }                                                                                                 \
    }                                                                                                     \
                                                                                                          \
    TEST_FIXTURE(PolyMultiplicationTests55, VEC##_##MUL##_mul_scal_div)                                   \
    {                                                                                                     \
        auto lhs = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>(1000000);             \
        auto rhs = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>(2000000);             \
                                                                                                          \
        rational_type scalar(2);                                                                          \
                                                                                                          \
        lhs.mul_scal_div(rhs, scalar);                                                                    \
        REQUIRE CHECK_EQUAL(3906, lhs.size());                                                            \
        for (auto item : lhs) {                                                                           \
            REQUIRE CHECK_EQUAL(construct_expected(item.key(), 1000000, 2000000) / 2, item.value());      \
        }                                                                                                 \
    }                                                                                                     \
                                                                                                          \
    TEST_FIXTURE(PolyMultiplicationTests55, VEC##_##MUL##_mul_scal_div_max_deg)                           \
    {                                                                                                     \
        auto lhs = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>(1000000);             \
        auto rhs = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>(2000000);             \
                                                                                                          \
        rational_type scalar(2);                                                                          \
                                                                                                          \
        lhs.mul_scal_div(rhs, scalar, DEG(3));                                                            \
        REQUIRE CHECK_EQUAL(3906, lhs.size());                                                            \
        for (auto item : lhs) {                                                                           \
            if (basis.degree(item.key()) <= 3) {                                                          \
                REQUIRE CHECK_EQUAL(construct_expected(item.key(), 1000000, 2000000) / 2, item.value());  \
            }                                                                                             \
            else {                                                                                        \
                REQUIRE CHECK_EQUAL(to_poly_key(item.key(), 1000000), item.value());                      \
            }                                                                                             \
        }                                                                                                 \
    }                                                                                                     \
                                                                                                          \
    TEST_FIXTURE(PolyMultiplicationTests55, VEC##_##MUL##_exp_test)                                       \
    {                                                                                                     \
        auto arg = generic_tensor<MUL##_multiplication, alg::vectors::VEC##_vector>();                    \
        arg[key_type()] = scalar_type();                                                                  \
                                                                                                          \
        tensor_type<MUL##_multiplication, alg::vectors::VEC##_vector> result{key_type(), scalar_type(1)}; \
        tensor_type<MUL##_multiplication, alg::vectors::VEC##_vector> one{key_type(), scalar_type(1)};    \
        for (DEG d = depth; d > 0; --d) {                                                                 \
            result.mul_scal_div(arg, rational_type(d), 5 - d + 1);                                        \
            result += one;                                                                                \
        }                                                                                                 \
                                                                                                          \
        alg::free_tensor<coeff_ring, width, depth, alg::vectors::VEC##_vector, alg::traditional_free_tensor_multiplication> expected_a(arg);                   \
        auto expected = exp(expected_a);                                                                  \
                                                                                                          \
        auto it = expected.cbegin();                                                                      \
        REQUIRE CHECK_EQUAL(3906, result.size());                                                          \
        REQUIRE CHECK_EQUAL(3906, expected.size());\
        for (auto item : result) {                                                                        \
            REQUIRE CHECK_EQUAL(it->key(), item.key());                                                   \
            REQUIRE CHECK_EQUAL(it->value().size(), item.value().size());                                 \
            if (it->value().size() != item.value().size()) {                                              \
                std::cout << it->key() << ' ' << item.value() - it->value() << "\n\n";                    \
            }                                                                                             \
            REQUIRE CHECK_EQUAL(it->value(), item.value());                                               \
            ++it;                                                                                         \
        }                                                                                                 \
    }                                                                                                     \
                                                                                                          \
    LA_TESTING_TENSOR_MAX_DEG_FMA(MUL, VEC, 0)                                                            \
    LA_TESTING_TENSOR_MAX_DEG_FMA(MUL, VEC, 1)                                                            \
    LA_TESTING_TENSOR_MAX_DEG_FMA(MUL, VEC, 2)                                                            \
    LA_TESTING_TENSOR_MAX_DEG_FMA(MUL, VEC, 3)                                                            \
    LA_TESTING_TENSOR_MAX_DEG_FMA(MUL, VEC, 4)                                                            \
    LA_TESTING_TENSOR_MAX_DEG_MUL_INPLACE(MUL, VEC, 0)                                                    \
    LA_TESTING_TENSOR_MAX_DEG_MUL_INPLACE(MUL, VEC, 1)                                                    \
    LA_TESTING_TENSOR_MAX_DEG_MUL_INPLACE(MUL, VEC, 2)                                                    \
    LA_TESTING_TENSOR_MAX_DEG_MUL_INPLACE(MUL, VEC, 3)                                                    \
    LA_TESTING_TENSOR_MAX_DEG_MUL_INPLACE(MUL, VEC, 4)                                                    \
    LA_TESTING_TENSOR_MAX_DEG_MSD(MUL, VEC, 0)                                                            \
    LA_TESTING_TENSOR_MAX_DEG_MSD(MUL, VEC, 1)                                                            \
    LA_TESTING_TENSOR_MAX_DEG_MSD(MUL, VEC, 2)                                                            \
    LA_TESTING_TENSOR_MAX_DEG_MSD(MUL, VEC, 3)                                                            \
    LA_TESTING_TENSOR_MAX_DEG_MSD(MUL, VEC, 4)                                                            \
    LA_TESTING_TENSOR_FMA_LOWER_DEG(MUL, VEC, 5, 3, 2)                                                    \
    LA_TESTING_TENSOR_FMA_LOWER_DEG(MUL, VEC, 5, 1, 4)                                                    \
    LA_TESTING_TENSOR_FMA_LOWER_DEG(MUL, VEC, 5, 4, 1)

    LA_TESTING_TENSOR_MUL_MAKE_TESTS(traditional, dense)
    LA_TESTING_TENSOR_MUL_MAKE_TESTS(tiled, dense)
    LA_TESTING_TENSOR_MUL_MAKE_TESTS(tiled2, dense)

    TEST_FIXTURE(PolyMultiplicationTests55, test_reverse_data_dense_extern)
    {
        auto lhs = generic_d_free_tensor(1000000);
        auto rhs = generic_d_free_tensor(2000000);

        const auto result = lhs * rhs;

        const auto& result_base = result.base_vector();
        const auto& reverse_data = result_base.reverse_data();

        CHECK(static_cast<bool>(reverse_data));
        for (auto item : result) {
            if (item.key().size() == depth) {
                break;
            }
            auto reverse_index = tensor_basis::key_to_index(item.key().reverse());
            REQUIRE CHECK_EQUAL(item.value(), reverse_data[reverse_index]);
        }
    }

    TEST_FIXTURE(PolyMultiplicationTests55, test_reverse_data_dense_extern_with_reverse)
    {
        auto lhs = generic_d_free_tensor(1000000);
        auto rhs = generic_d_free_tensor(2000000);

        lhs.base_vector().construct_reverse_data(depth - 1);

        const auto result = lhs * rhs;

        const auto& result_base = result.base_vector();
        const auto& reverse_data = result_base.reverse_data();

        CHECK(static_cast<bool>(reverse_data));
        for (auto item : result) {
            if (item.key().size() == depth) {
                break;
            }
            auto reverse_index = tensor_basis::key_to_index(item.key().reverse());
            REQUIRE CHECK_EQUAL(item.value(), reverse_data[reverse_index]);
        }
    }

    TEST_FIXTURE(PolyMultiplicationTests55, test_reverse_data_dense_inplace)
    {
        auto lhs = generic_d_free_tensor(1000000);
        auto rhs = generic_d_free_tensor(2000000);

        lhs *= rhs;

        const auto& result_base = lhs.base_vector();
        const auto& reverse_data = result_base.reverse_data();

        CHECK(static_cast<bool>(reverse_data));
        for (auto item : lhs) {
            if (item.key().size() == depth) {
                break;
            }
            auto reverse_index = tensor_basis::key_to_index(item.key().reverse());
            CHECK_EQUAL(construct_expected(item.key(), 1000000, 2000000), item.value());
            REQUIRE CHECK_EQUAL(item.value(), reverse_data[reverse_index]);
        }
    }

    TEST_FIXTURE(PolyMultiplicationTests55, test_reverse_data_dense_inplace_with_reverse)
    {
        auto lhs = generic_d_free_tensor(1000000);
        auto rhs = generic_d_free_tensor(2000000);

        lhs.base_vector().construct_reverse_data(depth - 1);

        lhs *= rhs;

        const auto& result_base = lhs.base_vector();
        const auto& reverse_data = result_base.reverse_data();

        CHECK(static_cast<bool>(reverse_data));
        for (auto item : lhs) {
            if (item.key().size() == depth) {
                break;
            }
            auto reverse_index = tensor_basis::key_to_index(item.key().reverse());
            REQUIRE CHECK_EQUAL(item.value(), reverse_data[reverse_index]);
        }
    }

    using PolyMultiplicationTests1283 = PolyMultiplicationTests<128, 3>;

#if 0
    TEST_FIXTURE(PolyMultiplicationTests1283, test_dense_mul_large_width)
    {
        auto lhs = generic_d_free_tensor(1000000);
        auto rhs = generic_d_free_tensor(2000000);

        using mul_type = alg::tiled_free_tensor_multiplication<width, depth, -2>;

        mul_type mul;
        tensor_type<mul_type, alg::vectors::dense_vector> result;
        alg::dtl::multiplication_traits<mul_type>::
                multiply_and_add(mul, result, lhs, rhs);

        CHECK_EQUAL(1 + width * (1 + width * (1 + width)), result.size());
        for (auto item : result) {
            REQUIRE CHECK_EQUAL(construct_expected(item.key(), 1000000, 2000000), item.value());
        }
    }

    using PolyMultiplicationTests1253 = PolyMultiplicationTests<125, 3>;
    TEST_FIXTURE(PolyMultiplicationTests1253, test_dense_mul_large_width_unequal_tiles)
    {
        auto lhs = generic_d_free_tensor(1000000);
        auto rhs = generic_d_free_tensor(2000000);

        using mul_type = alg::tiled_free_tensor_multiplication<width, depth, -2>;

        mul_type mul;
        tensor_type<mul_type, alg::vectors::dense_vector> result;
        alg::dtl::multiplication_traits<mul_type>::
                multiply_and_add(mul, result, lhs, rhs);

        CHECK_EQUAL(1 + width * (1 + width * (1 + width)), result.size());
        //        for (auto item : result) {
        //            REQUIRE CHECK_EQUAL(construct_expected(item.key(), 1000000, 2000000), item.value());
        //        }
        for (auto key : basis.iterate_keys()) {
            REQUIRE CHECK_EQUAL(construct_expected(key, 1000000, 2000000), result[key]);
        }
    }
#endif
}
