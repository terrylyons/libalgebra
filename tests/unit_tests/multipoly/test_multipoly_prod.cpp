//
// Created by user on 11/05/22.
//

#include <UnitTest++.h>
#include <libalgebra/libalgebra.h>
#include <libalgebra/alg_types.h>

SUITE(multipoly) {

    struct Fixture : alg_types<5, 5, Rational>
    {

    };


    TEST_FIXTURE(Fixture, prod_single_letters) {
        MULTIPOLY p1(LET(1), MULTIPOLY::one);
        MULTIPOLY p2(LET(2), MULTIPOLY::one);

        auto p = p1*p2;

        MULTIPOLY expected(typename MULTIPOLY::KEY{LET(1), LET(2)});

        CHECK_EQUAL(expected, p);
    }




}
