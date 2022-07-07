//
// Created by user on 07/07/22.
//

#include "test_tangents.h"


using namespace alg;

SUITE(tangents) {

struct TangentsFixture
{
    static constexpr DEG WIDTH = 5;
    static constexpr DEG DEPTH = 5;
    using coeff_type = coefficients::rational_field;
    using scalar_type = typename coeff_type::S;

    using tensor_basis_type = free_tensor_basis<WIDTH, DEPTH>;
    using tensor_t = free_tensor<coeff_type, WIDTH, DEPTH>;
    using tangent_tensor_t = tangent_vector<tensor_t>;

    using coeff_dist_t = la_testing::uniform_rational_distribution<typename coeff_type::S>;

    std::mt19937 rng;
    coeff_dist_t rational_dist;
    la_testing::random_vector_generator<tensor_t, coeff_dist_t> rvg;

    TangentsFixture() : rational_dist(-1, 1), rvg(rational_dist)
    {
        std::random_device dev;
        rng = std::mt19937(dev());
    }

    tensor_t random_tensor()
    {
        return rvg(rng);
    }

    tangent_tensor_t random_tangent_tensor()
    {
        return {random_tensor(), random_tensor()};
    }



};



    TEST_FIXTURE(TangentsFixture, test_tangent_multiplication_both_tangent_tensor)
    {
        const auto tt1 = random_tangent_tensor();
        const auto tt2 = random_tangent_tensor();

        const auto& rt1 = static_cast<const tensor_t&>(tt1);
        const auto& rt2 = static_cast<const tensor_t&>(tt2);

        auto result = tt1*tt2;

        auto expected_tensor = rt1*rt2;
        auto expected_tangent = rt1*tt2.tangent() + tt1.tangent()*rt2;

        CHECK_EQUAL(expected_tensor, static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected_tangent, result.tangent());
    }

    TEST_FIXTURE(TangentsFixture, test_tangent_multiplication_right_tensor)
    {
        const auto tt = random_tangent_tensor();
        const auto t = random_tensor();

        const auto& rt = static_cast<const tensor_t&>(tt);

        auto result = tt * t;

        auto expected_tensor = rt * t;
        auto expected_tangent = tt.tangent() * t;

        CHECK_EQUAL(expected_tensor, static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected_tangent, result.tangent());
    }

    TEST_FIXTURE(TangentsFixture, test_tangent_multiplication_left_tensor)
    {
        const auto t = random_tensor();
        const auto tt = random_tangent_tensor();

        const auto& rt = static_cast<const tensor_t&>(tt);

        auto result = t * tt;

        auto expected_tensor = t * rt;
        auto expected_tangent = t * tt.tangent();

        CHECK_EQUAL(expected_tensor, static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected_tangent, result.tangent());
    }

    TEST_FIXTURE(TangentsFixture, test_inplace_multiplication_tangent_tensor)
    {
        auto tt1 = random_tangent_tensor();
        const auto tt2 = random_tangent_tensor();

        const auto& rt1 = static_cast<const tensor_t&>(tt1);
        const auto& rt2 = static_cast<const tensor_t&>(tt2);

        auto expected_tensor = rt1*rt2;
        auto expected_tangent = rt1*tt2.tangent() + tt1.tangent()*rt2;

        tt1 *= tt2;

        CHECK_EQUAL(expected_tensor, static_cast<const tensor_t&>(tt1));
//        CHECK_EQUAL(expected_tangent, tt1.tangent());
    }

    TEST_FIXTURE(TangentsFixture, test_inplace_multiplication_tensor)
    {
        auto tt = random_tangent_tensor();
        auto t = random_tensor();

        const auto& rt = static_cast<const tensor_t&>(tt);

        auto expected_tensor = rt * t;
        auto expected_tangent = tt.tangent()*t;

        tt *= t;

        CHECK_EQUAL(expected_tensor, static_cast<const tensor_t&>(tt));
        CHECK_EQUAL(expected_tangent, tt.tangent());
    }


    TEST_FIXTURE(TangentsFixture, test_exp_tangent)
    {
        const auto tt = random_tangent_tensor();

        const auto& rt = static_cast<const tensor_t&>(tt);
        const auto& rtan = tt.tangent();

        DEG k;
        auto ad_X_k = [&](const tensor_t & a) {
            tensor_t result(a);
            for (DEG i=0; i<k; ++i) {
                result = rt*result - result*rt;
            }
            return result;
        };

        auto expected_tensor = exp(rt);
        tensor_t expected_tangent(expected_tensor);
        {
            tensor_t tmp(rtan);
            scalar_type factorial(1);
            for (k=1; k<=DEPTH; ++k) {
                factorial /= scalar_type(k+1);
                if (k % 2 == 0) {
                    tmp.add_scal_div(ad_X_k(rtan), factorial);
                } else {
                    tmp.sub_scal_div(ad_X_k(rtan), factorial);
                }
            }
            expected_tangent *= tmp;
        }

        auto result = exp(tt);

        CHECK_EQUAL(expected_tensor, static_cast<const tensor_t&>(result));
        CHECK_EQUAL(expected_tangent, result.tangent());
    }



}
