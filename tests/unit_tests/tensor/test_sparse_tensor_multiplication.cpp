//
// Created by sam on 23/03/2021.
//



#include <vector>
#include <map>

#include <UnitTest++.h>

#include <libalgebra/libalgebra.h>
#include <libalgebra/alg_types.h>

#include "../../common/time_and_details.h"
#include "../../common/rng.h"

using alg::LET;
using alg::DEG;

SUITE (sparse_tensor_multiplication) {

struct Fixture : public alg_types<5, 5, Rational>
{
    typedef alg_types<5, 5, Rational> ALG_TYPES;

    typedef typename ALG_TYPES::TENSOR::BASIS TBASIS;
    typedef typename TBASIS::KEY KEY;

    typedef typename ALG_TYPES::COEFF field;

    typedef alg::vectors::sparse_vector<
            TBASIS,
            field,
            std::map<KEY, typename field::S>
    > SPARSE;
    typedef alg::algebra<TBASIS, field,
                         alg::free_tensor_multiplication<5, 5>,
            alg::vectors::sparse_vector> TENSOR;

    const TENSOR tunit;
    const TENSOR tzero;

    Fixture() : tunit(KEY()), tzero()
    {}

    KEY make_key(const LET* arg, const std::size_t N)
    {
        KEY k;
        for (std::size_t i=0; i<N; ++i) {
            k.push_back(arg[i]);
        }
        return k;
    }

};

#define ADD_KEY(N, ...) \
    {                     \
        LET tmp[N] = {__VA_ARGS__};  \
        expected += TENSOR(make_key(tmp, N));\
    }


TEST_FIXTURE(Fixture, test_product_tunit_zero) {
    TEST_DETAILS();
    TENSOR lhs(tunit), rhs, expected;

            CHECK_EQUAL(expected, lhs * rhs);
}

TEST_FIXTURE(Fixture, test_product_zero_tunit) {
    TEST_DETAILS();
    TENSOR lhs, rhs(tunit), expected;

            CHECK_EQUAL(expected, lhs * rhs);
}


TEST_FIXTURE(Fixture, test_product_tunit_tunit) {
    TEST_DETAILS();

    TENSOR lhs(tunit), rhs(tunit);

            CHECK_EQUAL(tunit, lhs * rhs);
}

TEST_FIXTURE(Fixture, test_product_unidim_deg_1_tunit) {
    TEST_DETAILS();

    LET k1[] = {1};

    TENSOR lhs(tunit), rhs(make_key(k1, 1));

    TENSOR expected(rhs);
            CHECK_EQUAL(expected, lhs * rhs);
}

TEST_FIXTURE(Fixture, test_product_unidim_deg_1_1) {
    TEST_DETAILS();

    LET k1[] = {1};
    LET k2[] = {2};
    TENSOR lhs(make_key(k1, 1)), rhs(make_key(k2, 1));

    LET ek[] = {1, 2};
    TENSOR expected(make_key(ek,2));

            CHECK_EQUAL(expected, lhs * rhs);
}

TEST_FIXTURE(Fixture, test_product_unidim_deg_1_deg_2) {
    TEST_DETAILS();

    LET k1[] = {1};
    LET k2[] = {2, 3};
    TENSOR lhs(make_key(k1, 1)), rhs(make_key(k2, 2));

    LET ek[] = {1, 2, 3};
    TENSOR expected(make_key(ek, 3));

            CHECK_EQUAL(expected, lhs * rhs);
}

TEST_FIXTURE(Fixture, test_product_unidim_deg_2_deg_1) {
    TEST_DETAILS();

    LET k1[] = {1, 2};
    LET k2[] = {3};
    TENSOR lhs(make_key(k1, 2)), rhs(make_key(k2, 1));

    LET ek[] = {1, 2, 3};
    TENSOR expected(make_key(ek, 3));

            CHECK_EQUAL(expected, lhs * rhs);
}


TEST_FIXTURE(Fixture, test_product_unidim_deg_3_deg_3_overflow_removed) {
    TEST_DETAILS();

    LET k1[] = {1, 2, 3};
    LET k2[] = {1, 2, 3};
    TENSOR lhs(make_key(k1, 3)), rhs(make_key(k2, 3));

    TENSOR expected(tzero);

            CHECK_EQUAL(expected, lhs * rhs);
}

TEST_FIXTURE(Fixture, test_product_multiple_terms) {
    TEST_DETAILS();

    TENSOR lhs(tunit), rhs(tunit);

    LET letters[] = {1, 2, 3, 4, 4};
    // 1() + 1(1) + 1(1,2) + 1(1,2,3) + 1(1,2,3,4)
    for (DEG d=1; d<5; ++d) {
        lhs += TENSOR(make_key(&letters[0], d));
    }

    // 1() + 1(4) + 1(3,4) + 1(2,3,4) + 1(1,2,3,4)
    for (DEG d=1; d<5; ++d) {
        rhs += TENSOR(make_key(&letters[4-d], d));
    }

    TENSOR expected(tunit);


    const S one (1);

    // DEG 5
    ADD_KEY(5, 1,2,3,4,4);
    ADD_KEY(5, 1,2,3,3,4);
    ADD_KEY(5, 1,2,2,3,4);
    ADD_KEY(5, 1,1,2,3,4);

    /*
    expected += TENSOR(make_key((LET[]) {1,2,3,4,4}, 5));
    expected += TENSOR(make_key((LET[]) {1,2,3,3,4}, 5));
    expected += TENSOR(make_key((LET[]) {1,2,2,3,4}, 5));
    expected += TENSOR(make_key((LET[]) {1,1,2,3,4}, 5));
    */
    // DEG 4
    // (1,2,3,4) = () * (1,2,3,4), (1) * (2,3,4), (1,2) * (3,4), (1,2,3) * (4), (1,2,3,4) * ()
    {
        LET tmp[] = {1,2,3,4};
        expected.add_scal_prod(make_key(tmp, 4), S(5));
    }

    // DEG 3
    ADD_KEY(3, 1,2,3);
    ADD_KEY(3, 1,2,4);
    ADD_KEY(3, 1,3,4);
    ADD_KEY(3, 2,3,4);
    /*
    expected += TENSOR(make_key((LET[]) {1,2,3}, 3));
    expected += TENSOR(make_key((LET[]) {1,2,4}, 3));
    expected += TENSOR(make_key((LET[]) {1,3,4}, 3));
    expected += TENSOR(make_key((LET[]) {2,3,4}, 3));
*/

    // DEG 2
    ADD_KEY(2, 1,2);
    ADD_KEY(2, 1,4);
    ADD_KEY(2, 3,4);
    /*
    expected += TENSOR(make_key((LET[]) {1,2}, 2));
    expected += TENSOR(make_key((LET[]) {1,4}, 2));
    expected += TENSOR(make_key((LET[]) {3,4}, 2));
*/

    // DEG 1
    expected += TENSOR(KEY(LET(1)));
    expected += TENSOR(KEY(LET(4)));

#undef ADD_KEY

            CHECK_EQUAL(expected, lhs * rhs);
}




}
