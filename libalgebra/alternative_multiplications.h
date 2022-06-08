//
// Created by sam on 08/12/2021.
//

#ifndef LIBALGEBRA_ALTERNATIVE_MULTIPLICATIONS_H
#define LIBALGEBRA_ALTERNATIVE_MULTIPLICATIONS_H

#include <libalgebra/vectors/vectors.h>
#include <libalgebra/tensor.h>
#include <libalgebra/multiplication_helpers.h>
#include "half_shuffle_tensor_multiplication.h"
#include "area_tensor_multiplication.h"

namespace alg {
namespace dtl {

template<
        typename Multiplication,
        typename Op,
        typename Basis,
        typename Coeffs,
        template<typename, typename, typename...> class Vector,
        typename... Args>
void multiply_into_impl(
        const Multiplication& mpl, Op op,
        alg::vectors::vector<Basis, Coeffs, Vector, Args...>& result,
        const alg::vectors::vector<Basis, Coeffs, Vector, Args...>& lhs,
        const alg::vectors::vector<Basis, Coeffs, Vector, Args...>& rhs)
{
    mpl.multiply_and_add(result, lhs, rhs, op);
}

template<
        typename Multiplication,
        typename Op,
        typename Basis,
        typename Coeffs,
        template<typename, typename, typename...> class ArgVector,
        typename... ArgArgs,
        template<typename, typename, typename...> class ResultVector,
        typename... ResultArgs>
void multiply_into_impl(
        const Multiplication& mpl, Op op,
        alg::vectors::vector<Basis, Coeffs, ResultVector, ResultArgs...>& result,
        const alg::vectors::vector<Basis, Coeffs, ArgVector, ArgArgs...>& lhs,
        const alg::vectors::vector<Basis, Coeffs, ArgVector, ArgArgs...>& rhs)
{
    alg::vectors::vector<Basis, Coeffs, ArgVector, ArgArgs...> tmp;
    mpl.multiply_and_add(tmp, lhs, rhs, op);

    //// This is slow, but at the moment it is the only way add different vector
    //// types.
    //result.add_scal_prod(r.key(), r.value());
}

template<typename Multiplication, typename Op, typename Arg, typename Result = Arg>
void multiply_into(Result& result, const Arg& lhs, const Arg& rhs, const Multiplication& mpl, Op op)
{
    multiply_into_impl(mpl, op, result, lhs, rhs);
}

}



template<typename Tensor>
Tensor free_multiply(const Tensor& lhs, const Tensor& rhs)
{
    Tensor result;
    alg::free_tensor_multiplication<typename Tensor::coefficient_field> mpl;
    dtl::multiply_into(result, lhs, rhs, mpl, alg::mult::scalar_passthrough());
    return result;
}

template<typename Tensor>
Tensor shuffle_multiply(const Tensor& lhs, const Tensor& rhs)
{
    Tensor result;
    alg::shuffle_tensor_multiplication<typename Tensor::coefficient_field> mpl;
    dtl::multiply_into(result, lhs, rhs, mpl, alg::mult::scalar_passthrough());
    return result;
}

template<typename Tensor>
Tensor half_shuffle_multiply(const Tensor& lhs, const Tensor& rhs)
{
    Tensor result;
    half_shuffle_tensor_multiplication<typename Tensor::coefficient_field> mpl;
    dtl::multiply_into(result, lhs, rhs, mpl, alg::mult::scalar_passthrough());
    return result;
}

template<typename Tensor>
Tensor area_multiply(const Tensor& lhs, const Tensor& rhs)
{
    Tensor result;
    area_tensor_multiplication<typename Tensor::coefficient_field> mpl;
    dtl::multiply_into(result, lhs, rhs, mpl, alg::mult::scalar_passthrough());
    return result;
}

} // namespace alg


#endif //LIBALGEBRA_ALTERNATIVE_MULTIPLICATIONS_H
