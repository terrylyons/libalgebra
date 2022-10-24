//
// Created by sam on 23/03/2021.
//

#include <vector>

#include <UnitTest++.h>

#include <libalgebra/libalgebra.h>
#include <libalgebra/alg_types.h>

#include "../../common/time_and_details.h"
#include "../../common/helpers.h"



template<typename Coeff, unsigned Width, unsigned Depth>
struct fixture {

    typedef typename Coeff::S S;
    typedef typename Coeff::Q Q;

    typedef alg::free_tensor_basis<Width, Depth> TBASIS;
    typedef alg::vectors::dense_vector<TBASIS, Coeff> VECT;
    typedef alg::free_tensor<Coeff, Width, Depth, alg::vectors::dense_vector> TENSOR;

    typedef typename TBASIS::KEY KEY;

    template<class T>
    T exp_to_depth(T x, S one)
    {
        T result(one);
        T tone(S(1));

        for (size_t i = Depth; i >= 1; --i) {
            result *= (x / Q(i));
            result += tone;
        }
        return result;
    }

    TENSOR log_to_depth(const TENSOR& x)
    {
        TENSOR result;

        KEY kunit;
        TENSOR tunit(kunit);

        TENSOR xx(x);
        typename TENSOR::iterator it(xx.find(KEY()));
        if (it != xx.end()) {
            xx.erase(it);
        }

        for (unsigned i = Depth; i >= 1; --i) {

            if (i % 2 == 0) {
                result -= (tunit / Q(i));
            }
            else {
                result += (tunit / Q(i));
            }
            result *= xx;
        }
        return result;
    }
};


SUITE(dense_tensor_functions_float) {

    typedef alg::coefficients::coefficient_field<float> float_field;

    static const double expected_error = 3e-8f;
    typedef fixture<float_field, 5, 5> fixture;

#include "tensor_exp_suite.h"
#include "tensor_log_suite.h"
#include "tensor_inverse_suite.h"
#include "tensor_roundtrip_suite.h"


}


SUITE(dense_tensor_functions_double) {

    typedef alg::coefficients::coefficient_field<double> double_field;

    typedef fixture<double_field, 5, 5> fixture;
    static const double expected_error = 2e-15;

#include "tensor_exp_suite.h"
#include "tensor_log_suite.h"
#include "tensor_inverse_suite.h"
#include "tensor_roundtrip_suite.h"
}




SUITE(dense_tensor_functions_rational) {

    typedef typename alg_types<2, 2, Rational>::SCA Scalar;
    typedef alg::coefficients::coefficient_field<Scalar> rational_field;

    typedef fixture<rational_field, 5, 5> fixture;
    static const double expected_error = 0.0;

#include "tensor_exp_suite.h"
#include "tensor_log_suite.h"
#include "tensor_inverse_suite.h"
#include "tensor_roundtrip_suite.h"
}
