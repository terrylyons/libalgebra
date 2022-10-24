//
// Created by user on 08/06/22.
//

#include <UnitTest++.h>
#include <numeric>
#include <set>
#include <type_traits>


#include <libalgebra/libalgebra.h>
#include "half_shuffle_fixture.h"


using namespace alg;

constexpr DEG WIDTH = 3;
constexpr DEG DEPTH = 4;

constexpr DEG D = DEPTH;
constexpr DEG W = WIDTH;

using FixtureA = HalfShuffleFixture<DEPTH, WIDTH, Rational>;

SUITE(test_shuffles)
{
    // templated functions exist to allow cross type calculations
    // alg::free_multiply
    // alg::shuffle_multiply
    // alg::area_multiply
    // alg::half_shuffle_multiply

    TEST_FIXTURE(FixtureA, testhalfshuffle)
    {
        // tools to generate a random element

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> unif_dist(1, 1 << 5);

//        // create random LIE object
//        auto make_random_lie = [&](size_t d = D) {
//            LIE random_lie;
//            for (auto i = lbasis.begin(); i != lbasis.end() && lbasis.degree(i) <= d; i = lbasis.nextkey(i)) {
//                random_lie[i] = S(1) / S(unif_dist(gen));
//            }
//            return random_lie;
//        };
//
//        // create random HALF_SHUFFLE_TENSOR object
//        auto make_random_half_shuffle_tensor = [&](size_t d = D) {
//            HALF_SHUFFLE_TENSOR random_half_shuffle_tensor;
//            for (auto i = hsbasis.begin(); i != hsbasis.end() && hsbasis.degree(i) <= d; i = hsbasis.nextkey(i)) {
//                random_half_shuffle_tensor[i] = S(1) / S(unif_dist(gen));
//            }
//            return random_half_shuffle_tensor;
//        };
//
//        // create random SHUFFLE_TENSOR object
//        auto make_random_shuffle_tensor = [&](size_t d = D) {
//            SHUFFLE_TENSOR random_shuffle_tensor;
//            for (auto i = sbasis.begin(); i != sbasis.end() && sbasis.degree(i) <= d; i = sbasis.nextkey(i)) {
//                random_shuffle_tensor[i] = S(1) / S(unif_dist(gen));
//            }
//            return random_shuffle_tensor;
//        };

        // create random TENSOR object
        auto make_random_tensor = [&](size_t d = D) {
            TENSOR random_tensor;
            for (auto i = basis.begin(); i != basis.end() && basis.degree(i) <= d; i = basis.nextkey(i)) {
                random_tensor[i] = S(1) / S(unif_dist(gen));
            }
            return random_tensor;
        };

        for (auto i = basis.begin(); i != basis.end() && basis.degree(i) <= D; i = basis.nextkey(i))
            for (auto j = basis.begin(); j != i; j = basis.nextkey(j)) {
                TENSOR ti(i), tj(j);
                CHECK_EQUAL(TENSOR(), diff_shuffle_half_shuffle(ti, tj));
            }

        // define convenient tensors for tests
        TENSOR one(S(1)), random = make_random_tensor(), random1 = make_random_tensor();
        TENSOR t1(KEY(LET(1))), t2(KEY(LET(2))), t3(KEY(LET(3))), t11(t1 * t1), t12(t1 * t2), t13(t1 * t3), t21(t2 * t1), t22(t2 * t2), t23(t2 * t3), t31(t3 * t1), t32(t3 * t2), t33(t3 * t3);

        // check bi-linearity of half shuffle and the nilpotence

        TENSOR result, result2, lhs = make_random_tensor(), rhs = make_random_tensor();
        TENSOR::iterator lhscit = lhs.begin(), rhscit = rhs.begin();
        for (lhscit = lhs.begin(); lhscit != lhs.end(); ++lhscit)
            for (rhscit = rhs.begin(); rhscit != rhs.end(); ++rhscit) {

                TENSOR lhstmp;
                lhstmp[lhscit->key()] += lhscit->value();
                TENSOR rhstmp;
                rhstmp[rhscit->key()] += rhscit->value();
                if (basis.degree(lhscit->key()) + basis.degree(rhscit->key()) > D) {
                    result2 += alg::shuffle_multiply(lhstmp, rhstmp);
                }
                else
                    result += alg::shuffle_multiply(lhstmp, rhstmp);
            }
        CHECK_EQUAL(TENSOR(), result2);
        CHECK_EQUAL(result, alg::shuffle_multiply(lhs, rhs));
    }
    // correct on keys + bilinear means always correct

    TEST_FIXTURE(FixtureA, shuffle_powers)
    {
        TENSOR t1(KEY(LET(1))), t2(KEY(LET(2))), t3(KEY(LET(3))), t11(t1 * t1), t12(t1 * t2), t13(t1 * t3), t21(t2 * t1), t22(t2 * t2), t23(t2 * t3), t31(t3 * t1), t32(t3 * t2), t33(t3 * t3), t111(t1 * t11);
//        std::cout << alg::shuffle_multiply(t1, t1) << '\n' << t11*2 << "\n\n";
        CHECK_EQUAL(alg::shuffle_multiply(t1, t1), t11 * 2);
        CHECK_EQUAL(t111 * 3, alg::shuffle_multiply(t1, t11));
        CHECK_EQUAL(t111 * 3, alg::shuffle_multiply(t11, t1));
    }

    TEST_FIXTURE(FixtureA, shuffle_powers2)
    {
        SHUFFLE_TENSOR st1(KEY(LET(1))), st2(KEY(LET(2))), st3(KEY(LET(3))), st11(st1 * st1), st12(st1 * st2), st13(st1 * st3), st21(st2 * st1), st22(st2 * st2), st23(st2 * st3), st31(st3 * st1), st32(st3 * st2), st33(st3 * st3), st111(st1 * st11);
        CHECK_EQUAL(st11, alg::shuffle_multiply(st1, st1));
        CHECK_EQUAL(st111, alg::shuffle_multiply(st1, st11));
        CHECK_EQUAL(st111, alg::shuffle_multiply(st11, st1));
    }

    TEST_FIXTURE(FixtureA, testarea)
    {
        // tools to generate a random element

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> unif_dist(1, 1 << 5);

//        // create random LIE object
//        auto make_random_lie = [&](size_t d = D) {
//            LIE random_lie;
//            for (auto i = lbasis.begin(); i != lbasis.end() && lbasis.degree(i) <= d; i = lbasis.nextkey(i)) {
//                random_lie[i] = S(1) / S(unif_dist(gen));
//            }
//            return random_lie;
//        };
//
//        // create random HALF_SHUFFLE_TENSOR object
//        auto make_random_half_shuffle_tensor = [&](size_t d = D) {
//            HALF_SHUFFLE_TENSOR random_half_shuffle_tensor;
//            for (auto i = hsbasis.begin(); i != hsbasis.end() && hsbasis.degree(i) <= d; i = hsbasis.nextkey(i)) {
//                random_half_shuffle_tensor[i] = S(1) / S(unif_dist(gen));
//            }
//            return random_half_shuffle_tensor;
//        };
//
//        // create random SHUFFLE_TENSOR object
//        auto make_random_shuffle_tensor = [&](size_t d = D) {
//            SHUFFLE_TENSOR random_shuffle_tensor;
//            for (auto i = sbasis.begin(); i != sbasis.end() && sbasis.degree(i) <= d; i = sbasis.nextkey(i)) {
//                random_shuffle_tensor[i] = S(1) / S(unif_dist(gen));
//            }
//            return random_shuffle_tensor;
//        };
//
        // create random TENSOR object
        auto make_random_tensor = [&](size_t d = D) {
            TENSOR random_tensor;
            for (auto i = basis.begin(); i != basis.end() && basis.degree(i) <= d; i = basis.nextkey(i)) {
                random_tensor[i] = S(1) / S(unif_dist(gen));
            }
            return random_tensor;
        };


        for (auto i = basis.begin(); i != basis.end() && basis.degree(i) <= D; i = basis.nextkey(i))
            for (auto j = basis.begin(); j != i; j = basis.nextkey(j)) {
                TENSOR ti(i), tj(j);
                CHECK_EQUAL(TENSOR(), diff_area_shuffle(ti, tj));
            }

        // define convenient tensors for tests
        TENSOR one(S(1)), random = make_random_tensor(), random1 = make_random_tensor();
        TENSOR t1(KEY(LET(1))), t2(KEY(LET(2))), t3(KEY(LET(3))), t11(t1 * t1), t12(t1 * t2), t13(t1 * t3), t21(t2 * t1), t22(t2 * t2), t23(t2 * t3), t31(t3 * t1), t32(t3 * t2), t33(t3 * t3);

        // check bi-linearity of area and the nilpotence

        TENSOR result, result2, lhs = make_random_tensor(), rhs = make_random_tensor();
        TENSOR::iterator lhscit = lhs.begin(), rhscit = rhs.begin();
        for (lhscit = lhs.begin(); lhscit != lhs.end(); ++lhscit)
            for (rhscit = rhs.begin(); rhscit != rhs.end(); ++rhscit) {

                TENSOR lhstmp;
                lhstmp[lhscit->key()] += lhscit->value();
                TENSOR rhstmp;
                rhstmp[rhscit->key()] += rhscit->value();
                if (basis.degree(lhscit->key()) + basis.degree(rhscit->key()) > D) {
                    result2 += alg::area_multiply(lhstmp, rhstmp);
                }
                else
                    result += alg::area_multiply(lhstmp, rhstmp);
            }
        CHECK_EQUAL(TENSOR(), result2);
        CHECK_EQUAL(result, alg::area_multiply(lhs, rhs));
        // correct on keys + bilinear means always correct
    }
}
