//
// Created by sam on 16/11/2021.
//

#ifndef LIBALGEBRA_FUNCTIONALS_H
#define LIBALGEBRA_FUNCTIONALS_H

#include "implementation_types.h"

#include "dot_product_implementations.h"
#include "operators.h"
#include "tensor.h"
#include "scalar_bundle.h"
#include "vector_bundle.h"

namespace alg {
namespace operators {

template <typename Impl>
class linear_functional : public linear_operator<Impl> {
public:

    using linear_operator<Impl>::linear_operator;

    using Impl::operator();


    template <typename ArgBase, typename ArgFibre>
    scalar_bundle<typename ArgBase::coefficient_field, typename ArgFibre::coefficient_field>
    operator()(const vector_bundle<ArgBase, ArgFibre>& arg) const {
        return { (*this)(arg.base()), (*this)(arg.fibre()) };
    }


};


template<typename ShuffleTensor, typename FreeTensor>
using shuffle_tensor_functional = linear_functional<dot_product_implementation<ShuffleTensor, FreeTensor>>;


template <typename ShuffleTensor, typename FreeTensor>
typename FreeTensor::SCALAR apply_functional(const ShuffleTensor& functional, const FreeTensor& tensor)
{
    shuffle_tensor_functional<ShuffleTensor, FreeTensor> shfunctional(functional);
    return shfunctional(tensor);
}


}// namespace operators
}// namespace alg

#endif//LIBALGEBRA_FUNCTIONALS_H
