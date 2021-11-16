//
// Created by sam on 16/11/2021.
//

#ifndef LIBALGEBRA_FUNCTIONALS_H
#define LIBALGEBRA_FUNCTIONALS_H

#include <libalgebra/implementation_types.h>

#include <libalgebra/operators/dot_product_implementations.h>
#include <libalgebra/operators/operators.h>

#include <libalgebra/tensor.h>

namespace alg {
namespace operators {

template<typename ShuffleTensor, typename FreeTensor>
using shuffle_tensor_functional = linear_functional<dot_product_implementation<ShuffleTensor, FreeTensor>, FreeTensor>;

}// namespace operators
}// namespace alg

#endif//LIBALGEBRA_FUNCTIONALS_H
