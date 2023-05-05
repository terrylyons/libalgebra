//
// Created by user on 27/04/23.
//

#include <libalgebra/reconfigure_operator.h>
#include <libalgebra/rational_coefficients.h>


#include <UnitTest++/UnitTest++.h>

using namespace alg;

SUITE(recalibrate_operator) {

    TEST(TestRecalibrateFreeTensorDownDegree) {

        operators::recalibrate_operator<operators::change_depth<3>> op;

        free_tensor<alg::coefficients::rational_field, 5, 5, vectors::dense_vector, traditional_free_tensor_multiplication> src;
        for (DIMN i=0; i < src.basis.start_of_degree(6); ++i) {
            src.add_scal_prod(src.basis.index_to_key(i), coefficients::rational(i));
        }

        auto recal = op(src);

        CHECK_EQUAL(3, recal.degree());

    }

    TEST(TestRecalibrateFreeTensorUpDegree) {

        operators::recalibrate_operator<operators::change_depth<6>> op;

        free_tensor<alg::coefficients::rational_field, 5, 5, vectors::dense_vector, traditional_free_tensor_multiplication> src;
        for (DIMN i=0; i < src.basis.start_of_degree(6); ++i) {
            src.add_scal_prod(src.basis.index_to_key(i), coefficients::rational(i));
        }

        auto recal = op(src);

        CHECK_EQUAL(5, recal.degree());

    }



}
