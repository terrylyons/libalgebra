//
// Created by sam on 06/10/2021.
//

#include <libalgebra/libalgebra.h>
#include <UnitTest++/UnitTest++.h>
#include <libalgebra/alg_types.h>


using fixture55double = alg_types<5, 5, DPReal>;

SUITE(tensor_dual_pairing) {


    TEST_FIXTURE(fixture55double, test_bracket_zeros) {
        TENSOR ftensor;
        SHUFFLE_TENSOR stensor;

        CHECK_EQUAL(0.0, alg::apply_functional(stensor, ftensor));

    }





}