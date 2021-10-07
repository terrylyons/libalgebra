//
// Created by sam on 06/10/2021.
//

#include <vector>

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

    TEST_FIXTURE(fixture55double, test_bracket_identity) {
        TENSOR ftensor(S(1));
        SHUFFLE_TENSOR stensor(S(1));

        CHECK_EQUAL(1.0, alg::apply_functional(stensor, ftensor));
    }

    TEST_FIXTURE(fixture55double, test_bracket_deg1) {
        std::vector<double> fdata { 0.0, 1.0, 2.0, 3.0 };
        std::vector<double> sdata { 0.0, 1.0, 1.0, 1.0 };

        TENSOR ftensor(&*fdata.begin(), &*fdata.end());
        SHUFFLE_TENSOR stensor(&*sdata.begin(), &*sdata.end());

        CHECK_EQUAL(6.0, alg::apply_functional(stensor, ftensor));
    }

    TEST_FIXTURE(fixture55double, test_bracket_deg1_mismatch) {
        std::vector<double> fdata { 0.0, 1.0, 2.0, 3.0 };
        std::vector<double> sdata { 0.0, 0.0, 1.0, 0.0 };

        TENSOR ftensor(&*fdata.begin(), &*fdata.end());
        SHUFFLE_TENSOR stensor(&*sdata.begin(), &*sdata.end());

        CHECK_EQUAL(2.0, alg::apply_functional(stensor, ftensor));
    }

    TEST_FIXTURE(fixture55double, test_bracket_deg1_complete_mismatch) {
        std::vector<double> fdata { 0.0, 1.0, 2.0, 3.0 };
        std::vector<double> sdata { 1.0, 0.0, 0.0, 0.0 };

        TENSOR ftensor(&*fdata.begin(), &*fdata.end());
        SHUFFLE_TENSOR stensor(&*sdata.begin(), &*sdata.end());

        CHECK_EQUAL(0.0, alg::apply_functional(stensor, ftensor));
    }





}