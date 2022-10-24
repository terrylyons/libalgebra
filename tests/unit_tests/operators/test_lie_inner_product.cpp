//
// Created by sam on 31/01/2022.
//

#include <libalgebra/alg_types.h>
#include <libalgebra/libalgebra.h>
#include <libalgebra/lie_inner_product.h>

#include <UnitTest++.h>

SUITE(lie_inner_product_tests)
{

    struct fixture : alg_types<5, 5, Rational> {
        using LIE_IP = alg::operators::lie_inner_product<COEFF, ALPHABET_SIZE, DEPTH, alg::vectors::sparse_vector>;

        LIE_IP ip;
    };

    TEST_FIXTURE(fixture, test_zero_ip)
    {
        LIE left, right;

        CHECK_EQUAL(S(0), ip(left, right));
    }

    TEST_FIXTURE(fixture, test_ip_paired_letters)
    {

        for (LET l = 1; l <= ALPHABET_SIZE; ++l) {
            LIE left(l), right(l);

            CHECK_EQUAL(S(1), ip(left, right));
        }
    }

    TEST_FIXTURE(fixture, test_ip_different_letters)
    {
        for (LET l1 = 1; l1 <= ALPHABET_SIZE; ++l1) {
            for (LET l2 = 1; l2 <= ALPHABET_SIZE; ++l2) {
                if (l1 != l2) {
                    LIE left(l1), right(l2);
                    CHECK_EQUAL(S(0), ip(left, right));
                }
            }
        }
    }

    TEST_FIXTURE(fixture, test_ip_deg2_key)
    {
        typename LIE::KEY k1(21);// [2, [1, 3]]
        typename LIE::KEY k2(27);// [3, [1, 2]]

        LIE left(k1, S(1)), right(k2, S(1));
        S expected = 2;// from (2 tensor 1 tensor 3 and 3 tensor 1 tensor 2)

        std::cout << left << ' ' << right << '\n';
        CHECK_EQUAL(expected, ip(left, right));
    }
}
