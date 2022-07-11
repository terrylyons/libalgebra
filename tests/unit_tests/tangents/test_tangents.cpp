//
// Created by user on 07/07/22.
//

#include "test_tangents.h"


using namespace alg;

SUITE(tangents) {

struct TangentsFixture
{
    static constexpr DEG WIDTH = 5;
    static constexpr DEG DEPTH = 2;

    using rational_field = coefficients::rational_field;
    using rational = typename rational_field::S;

    using coeff_type = coefficients::rational_field;
//    using coeff_type = coefficients::coefficient_ring<poly<rational_field>, typename rational_field::Q>;

    using scalar_type = typename coeff_type::S;

    using tensor_basis_type = free_tensor_basis<WIDTH, DEPTH>;
    using tensor_t = free_tensor<coeff_type, WIDTH, DEPTH>;
    using tangent_tensor_t = vector_bundle<tensor_t>;

    using coeff_dist_t = la_testing::uniform_rational_distribution<typename rational_field::S>;

    std::mt19937 rng;
    coeff_dist_t rational_dist;
    la_testing::random_vector_generator<tensor_t, coeff_dist_t> rvg;

    TangentsFixture() : rational_dist(-1, 1), rvg(rational_dist)
    {
        std::random_device dev;
        rng = std::mt19937(dev());
    }

    tensor_t random_tensor(LET i=0)
    {
        return rvg(rng);
    }

//    tensor_t random_tensor(LET i=0)
//    {
//
//        tensor_t result;
//        LET ki=0;
//        tensor_basis_type basis;
//        for (auto k : basis.iterate_keys()) {
//            result.add_scal_prod(k, scalar_type(i + (ki++), rational(1)));
//        }
//        return result;
//    }

    tangent_tensor_t random_tangent_tensor()
    {
        return {random_tensor(100), random_tensor(200)};
    }

    scalar_type random_scalar(LET offset=1000)
    {
//        return scalar_type(offset, rational(1));
        return rational_dist(rng);
    }



};

    TEST_FIXTURE(TangentsFixture, test_equality_equals) {
        const auto tt = random_tangent_tensor();

        CHECK_EQUAL(tt, tt);
    }

    TEST_FIXTURE(TangentsFixture, test_equality_fails) {
        const auto tt1 = random_tangent_tensor();
        const auto tt2 = random_tangent_tensor();

        CHECK(!(tt1 == tt2));
    }

    TEST_FIXTURE(TangentsFixture, test_equality_fails_different_fibres) {
        const auto t = random_tensor();
        const auto fibre1 = random_tensor();
        const auto fibre2 = random_tensor();

        const tangent_tensor_t tt1(t, fibre1);
        const tangent_tensor_t tt2(t, fibre2);

        CHECK(!(tt1 == tt2));
    }

    TEST_FIXTURE(TangentsFixture, test_equality_fails_different_vectors) {
        const auto t1 = random_tensor();
        const auto t2 = random_tensor();
        const auto fibre = random_tensor();

        const tangent_tensor_t tt1(t1, fibre), tt2(t2, fibre);

        CHECK(!(tt1 == tt2));
    }

    TEST_FIXTURE(TangentsFixture, test_notequality_fails) {
        const auto tt = random_tangent_tensor();

        CHECK(!(tt != tt));
    }

    TEST_FIXTURE(TangentsFixture, test_notequality_true) {
        const auto tt1 = random_tangent_tensor();
        const auto tt2 = random_tangent_tensor();

        CHECK(tt1 != tt2);
    }

    TEST_FIXTURE(TangentsFixture, test_notequality_different_fibres) {
        const auto t = random_tensor();
        const auto fibre1 = random_tensor();
        const auto fibre2 = random_tensor();

        const tangent_tensor_t tt1(t, fibre1);
        const tangent_tensor_t tt2(t, fibre2);

        CHECK((tt1 != tt2));
    }

    TEST_FIXTURE(TangentsFixture, test_notequality_different_vectors) {
        const auto t1 = random_tensor();
        const auto t2 = random_tensor();
        const auto fibre = random_tensor();

        const tangent_tensor_t tt1(t1, fibre), tt2(t2, fibre);

        CHECK((tt1 != tt2));
    }


    TEST_FIXTURE(TangentsFixture, test_tangent_unary_minus)
    {
        const auto tt = random_tangent_tensor();

        auto result = -tt;
        auto expected_tensor = -static_cast<const tensor_t&>(tt);
        auto expected_tangent = -tt.fibre();

        CHECK_EQUAL(expected_tensor, static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected_tangent, result.fibre());
    }

    TEST_FIXTURE(TangentsFixture, test_addition_both_tangents)
    {
        const auto tt1 = random_tangent_tensor();
        const auto tt2 = random_tangent_tensor();

        auto result = tt1 + tt2;
        auto expected_tensor = static_cast<const tensor_t&>(tt1) + static_cast<const tensor_t&>(tt2);
        auto expected_tangent = tt1.fibre() + tt2.fibre();

        CHECK_EQUAL(expected_tensor, static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected_tangent, result.fibre());
    }

    TEST_FIXTURE(TangentsFixture, test_tangent_add_tensor_right)
    {
        const auto tt = random_tangent_tensor();
        const auto t = random_tensor();

        auto result = tt + t;
        auto expected_tensor = static_cast<const tensor_t&>(tt) + t;
        auto expected_tangent = tt.fibre();

        CHECK_EQUAL(expected_tensor, static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected_tangent, result.fibre());
    }

    TEST_FIXTURE(TangentsFixture, test_tangent_add_tensor_left)
    {
        const auto t = random_tensor();
        const auto tt = random_tangent_tensor();

        auto result = t + tt;
        auto expected_tensor = t + static_cast<const tensor_t&>(tt);
        auto expected_tangent = tt.fibre();

        CHECK_EQUAL(expected_tensor, static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected_tangent, result.fibre());
    }

    TEST_FIXTURE(TangentsFixture, test_tangent_subtraction)
    {
        const auto tt1 = random_tangent_tensor();
        const auto tt2 = random_tangent_tensor();

        auto result = tt1 - tt2;
        auto expected_tensor = static_cast<const tensor_t&>(tt1) - static_cast<const tensor_t&>(tt2);
        auto expected_tangent = tt1.fibre() - tt2.fibre();

        CHECK_EQUAL(expected_tensor, static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected_tangent, result.fibre());
    }

    TEST_FIXTURE(TangentsFixture, test_tangent_subtraction_right_tensor)
    {
        const auto tt = random_tangent_tensor();
        const auto t = random_tensor();

        auto result = tt + t;
        auto expected_tensor = static_cast<const tensor_t&>(tt) + t;
        auto expected_tangent = tt.fibre();

        CHECK_EQUAL(expected_tensor, static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected_tangent, result.fibre());
    }

TEST_FIXTURE(TangentsFixture, test_tangent_multiplication_both_tangent_tensor)
    {
        const auto tt1 = random_tangent_tensor();
        const auto tt2 = random_tangent_tensor();

        const auto& rt1 = static_cast<const tensor_t&>(tt1);
        const auto& rt2 = static_cast<const tensor_t&>(tt2);

        auto result = tt1*tt2;

        auto expected_tensor = rt1*rt2;
        auto expected_tangent = rt1* tt2.fibre() + tt1.fibre()*rt2;

        CHECK_EQUAL(expected_tensor, static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected_tangent, result.fibre());
    }

    TEST_FIXTURE(TangentsFixture, test_tangent_multiplication_right_tensor)
    {
        const auto tt = random_tangent_tensor();
        const auto t = random_tensor();

        const auto& rt = static_cast<const tensor_t&>(tt);

        auto result = tt * t;

        auto expected_tensor = rt * t;
        auto expected_tangent = tt.fibre() * t;

        CHECK_EQUAL(expected_tensor, static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected_tangent, result.fibre());
    }

    TEST_FIXTURE(TangentsFixture, test_tangent_multiplication_left_tensor)
    {
        const auto t = random_tensor();
        const auto tt = random_tangent_tensor();

        const auto& rt = static_cast<const tensor_t&>(tt);

        auto result = t * tt;

        auto expected_tensor = t * rt;
        auto expected_tangent = t * tt.fibre();

        CHECK_EQUAL(expected_tensor, static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected_tangent, result.fibre());
    }

    TEST_FIXTURE(TangentsFixture, test_inplace_multiplication_tangent_tensor)
    {
        auto tt1 = random_tangent_tensor();
        const auto tt2 = random_tangent_tensor();

        const auto& rt1 = static_cast<const tensor_t&>(tt1);
        const auto& rt2 = static_cast<const tensor_t&>(tt2);

        auto expected_tensor = rt1*rt2;
        auto expected_tangent = rt1* tt2.fibre() + tt1.fibre()*rt2;

        tt1 *= tt2;

        CHECK_EQUAL(expected_tensor, static_cast<const tensor_t&>(tt1));
        CHECK_EQUAL(expected_tangent, tt1.fibre());
    }

    TEST_FIXTURE(TangentsFixture, test_inplace_multiplication_tensor)
    {
        auto tt = random_tangent_tensor();
        auto t = random_tensor();

        const auto& rt = static_cast<const tensor_t&>(tt);

        auto expected_tensor = rt * t;
        auto expected_tangent = tt.fibre()*t;

        tt *= t;

        CHECK_EQUAL(expected_tensor, static_cast<const tensor_t&>(tt));
        CHECK_EQUAL(expected_tangent, tt.fibre());
    }

    TEST_FIXTURE(TangentsFixture, test_add_scal_prod_tangent)
    {
        auto tt1 = random_tangent_tensor();
        const auto tt2 = random_tangent_tensor();
        const auto rs = random_scalar();

        auto expected = tt1 + (tt2 * rs);
        auto result = tt1.add_scal_prod(tt2, rs);

        CHECK_EQUAL(static_cast<const tensor_t&>(expected), static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected.fibre(), result.fibre());
    }


    TEST_FIXTURE(TangentsFixture, test_add_scal_prod_tensor)
    {
        auto tt1 = random_tangent_tensor();
        const auto t2 = random_tensor();
        const auto rs = random_scalar();

        auto expected = tt1 + (t2 * rs);
        auto result = tt1.add_scal_prod(t2, rs);

        CHECK_EQUAL(static_cast<const tensor_t&>(expected), static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected.fibre(), result.fibre());
    }


    TEST_FIXTURE(TangentsFixture, test_sub_scal_prod_tangent)
    {
        auto tt1 = random_tangent_tensor();
        const auto tt2 = random_tangent_tensor();
        const auto rs = random_scalar();

        auto expected = tt1 - (tt2 * rs);
        auto result = tt1.sub_scal_prod(tt2, rs);

        CHECK_EQUAL(static_cast<const tensor_t&>(expected), static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected.fibre(), result.fibre());
    }


    TEST_FIXTURE(TangentsFixture, test_sub_scal_prod_tensor)
    {
        auto tt1 = random_tangent_tensor();
        const auto t2 = random_tensor();
        const auto rs = random_scalar();

        auto expected = tt1 - (t2 * rs);
        auto result = tt1.sub_scal_prod(t2, rs);

        CHECK_EQUAL(static_cast<const tensor_t&>(expected), static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected.fibre(), result.fibre());
    }

    TEST_FIXTURE(TangentsFixture, test_add_scal_div_tangent)
    {
        auto tt1 = random_tangent_tensor();
        const auto tt2 = random_tangent_tensor();
        const auto rs = rational_dist(rng);

        auto expected = tt1 + (tt2 / rs);
        auto result = tt1.add_scal_div(tt2, rs);

        CHECK_EQUAL(static_cast<const tensor_t&>(expected), static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected.fibre(), result.fibre());
    }


    TEST_FIXTURE(TangentsFixture, test_add_scal_div_tensor)
    {
        auto tt1 = random_tangent_tensor();
        const auto t2 = random_tensor();
        const auto rs = rational_dist(rng);

        auto expected = tt1 + (t2 / rs);
        auto result = tt1.add_scal_div(t2, rs);

        CHECK_EQUAL(static_cast<const tensor_t&>(expected), static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected.fibre(), result.fibre());
    }

    TEST_FIXTURE(TangentsFixture, test_sub_scal_div_tangent)
    {
        auto tt1 = random_tangent_tensor();
        const auto tt2 = random_tangent_tensor();
        const auto rs = rational_dist(rng);

        auto expected = tt1 - (tt2 / rs);
        auto result = tt1.sub_scal_div(tt2, rs);

        CHECK_EQUAL(static_cast<const tensor_t&>(expected), static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected.fibre(), result.fibre());
    }


    TEST_FIXTURE(TangentsFixture, test_sub_scal_div_tensor)
    {
        auto tt1 = random_tangent_tensor();
        const auto t2 = random_tensor();
        const auto rs = rational_dist(rng);

        auto expected = tt1 - (t2 / rs);
        auto result = tt1.sub_scal_div(t2, rs);

        CHECK_EQUAL(static_cast<const tensor_t&>(expected), static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected.fibre(), result.fibre());
    }


    TEST_FIXTURE(TangentsFixture, test_mul_scal_prod_tangent) {
        auto tt1 = random_tangent_tensor();
        const auto tt2 = random_tangent_tensor();
        const auto rs = random_scalar();

        auto expected = tt1 * (tt2 * rs);
        auto result = tt1.mul_scal_prod(tt2, rs);

        CHECK_EQUAL(static_cast<const tensor_t&>(expected), static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected.fibre(), result.fibre());
    }

    TEST_FIXTURE(TangentsFixture, test_mul_scal_prod_tensor) {
        auto tt1 = random_tangent_tensor();
        const auto t2 = random_tensor();
        const auto rs = random_scalar();

        auto expected = tt1 * (t2 * rs);
        auto result = tt1.mul_scal_prod(t2, rs);

        CHECK_EQUAL(static_cast<const tensor_t&>(expected), static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected.fibre(), result.fibre());
    }

    TEST_FIXTURE(TangentsFixture, test_mul_scal_div_tangent) {
        auto tt1 = random_tangent_tensor();
        const auto tt2 = random_tangent_tensor();
        const auto rs = rational_dist(rng);

        auto expected = tt1 * (tt2 / rs);
        auto result = tt1.mul_scal_div(tt2, rs);

        CHECK_EQUAL(static_cast<const tensor_t&>(expected), static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected.fibre(), result.fibre());
    }

    TEST_FIXTURE(TangentsFixture, test_mul_scal_div_tensor) {
        auto tt1 = random_tangent_tensor();
        const auto t2 = random_tensor();
        const auto rs = rational_dist(rng);

        auto expected = tt1 * (t2 / rs);
        auto result = tt1.mul_scal_div(t2, rs);

        CHECK_EQUAL(static_cast<const tensor_t&>(expected), static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected.fibre(), result.fibre());
    }


    TEST_FIXTURE(TangentsFixture, test_exp_tangent)
    {
        const auto tt = random_tangent_tensor();

        const auto& rt = static_cast<const tensor_t&>(tt);
        const auto& rtan = tt.fibre();

        typename tensor_t::KEY kunit;
        tensor_t tunit(kunit), expected_tensor(kunit), expected_tangent(kunit);
        for (DEG i = DEPTH; i >= 1; --i) {
            typename tensor_t::RATIONAL divisor(i);

            tensor_t& lhs = expected_tensor;
            tensor_t& lhs_tan = expected_tangent;


            expected_tangent = expected_tensor * (rtan / divisor) + expected_tangent * (rt / divisor);

            expected_tensor.mul_scal_div(rt, divisor);

            expected_tensor += tunit;
            expected_tangent += tunit;
        }

        auto result = exp(tt);

        CHECK_EQUAL(expected_tensor, static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected_tangent, result.fibre());
    }
//
//    TEST_FIXTURE(TangentsFixture, test_exp_via_formula) {
//
//        const auto tt = random_tangent_tensor();
//
//        const auto& rt = static_cast<const tensor_t&>(tt);
//        const auto& rtan = tt.fibre();
//
//        auto ad_X_k = [&](DEG k) {
//            tensor_t result(rtan);
//            for (DEG i = 0; i < k; ++i) {
//                result = rt * result - result * rt;
//            }
//            return result;
//        };
//
//        auto expected_tensor = exp(rt);
//        auto expected_tangent = expected_tensor;
//        {
//            tensor_t tmp(tt.fibre());
//            typename tensor_t::RATIONAL divisor(1);
//            for (DEG k=1; k<=DEPTH; ++k) {
//                divisor *= k+1;
//                if (k % 2 == 0) {
//                    tmp.add_scal_div(ad_X_k(k), divisor);
//                } else {
//                    tmp.sub_scal_div(ad_X_k(k), divisor);
//                }
//            }
//            expected_tangent *= tmp;
//        }
//
//    auto result = exp(tt);
//
//    CHECK_EQUAL(expected_tensor, static_cast<const tensor_t&>(result));
//    CHECK_EQUAL(expected_tangent, result.fibre());
//}


TEST_FIXTURE(TangentsFixture, tensor_exp_sanity_test) {

    const auto tt = random_tangent_tensor();

    const auto& rt = static_cast<const tensor_t&>(tt);
    const auto& rtan = tt.fibre();

    typename tensor_t::KEY kunit;
    tangent_tensor_t tunit(kunit);
    auto expected = tunit;
    for (DEG i=DEPTH; i>=1; --i) {
        expected.mul_scal_div(tt, typename tensor_t::RATIONAL(i));
        expected += tunit;
    }

    auto result = exp(tt);

    CHECK_EQUAL(expected, result);

}

}
