//
// Created by sam on 08/09/2021.
//

#include <vector>

#include <UnitTest++.h>

#include <libalgebra/libalgebra.h>
#include <libalgebra/coefficients.h>

#include "../../common/helpers.h"



SUITE(tensor_signature_calc_check) {

struct TensorSignatureCalcFixture {

    typedef alg::free_tensor_basis<2, 2> tensor_basis;
    typedef typename tensor_basis::KEY tkey;
    typedef alg::coefficients::double_field coeff_t;
    typedef alg::free_tensor_multiplication<2, 2> tensor_multiplication_t;

    typedef alg::vectors::dense_vector<tensor_basis, coeff_t> dense_vec_t;

    typedef alg::free_tensor<coeff_t, 2, 2, alg::vectors::dense_vector> dense_tensor;

    typedef alg::lie<coeff_t, 2, 1> lie_increment;
    typedef alg::maps<coeff_t, 2, 2, dense_tensor> maps;
    typedef alg::cbh<coeff_t, 2, 2, dense_tensor> cbh;

    std::vector<lie_increment> increments;

    dense_tensor expected;

    TensorSignatureCalcFixture(double lower = 0.0, double upper = 2.0, unsigned length = 5) : increments(), expected()
    {
        increments.reserve(length);

        lie_increment tmp;
        tmp.add_scal_prod(alg::LET(1), 2.0 / length);
        tmp.add_scal_prod(alg::LET(2), 4.0 / length);

        for (unsigned i=0; i<length; ++i) {
            increments.push_back(tmp);
        }

        expected.add_scal_prod(tkey(), 1.0);

        tkey k1(alg::LET(1)), k2(alg::LET(2));

        expected.add_scal_prod(k1, 2.0);
        expected.add_scal_prod(k2, 4.0);
        expected.add_scal_prod(k1*k1, 2.0);
        expected.add_scal_prod(k1*k2, 4.0);
        expected.add_scal_prod(k2*k1, 4.0);
        expected.add_scal_prod(k2*k2, 8.0);
    }

};

    TEST_FIXTURE(TensorSignatureCalcFixture, check_signature_calculation_correct_exp) {

        dense_tensor tmp(1);

        typename std::vector<lie_increment>::iterator it;
        maps maps_;

        for (it=increments.begin(); it != increments.end(); ++it) {
            tmp *= exp(maps_.l2t(*it));
        }

        CHECK_VEC_CLOSE(expected, tmp, 2e-15);

    }

    TEST_FIXTURE(TensorSignatureCalcFixture, check_signature_calculation_correct_fmexp) {

        dense_tensor tmp(1);

        typename std::vector<lie_increment>::iterator it;
        maps maps_;

        for (it=increments.begin(); it != increments.end(); ++it) {
            tmp.fmexp_inplace(maps_.l2t(*it));
        }

        CHECK_VEC_CLOSE(expected, tmp, 2e-15);
    }

    TEST_FIXTURE(TensorSignatureCalcFixture, check_log_exp_roundtrip) {

        CHECK_VEC_CLOSE(expected, exp(log(expected)), 2e-15);

    }

    TEST_FIXTURE(TensorSignatureCalcFixture, check_exp_cbh) {

        maps maps_;
        cbh cbh_;
        dense_tensor result = exp(maps_.l2t(cbh_.full(increments.begin(), increments.end())));

        CHECK_VEC_CLOSE(expected, result, 2e-15);

    }


}
