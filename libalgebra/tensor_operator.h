//
// Created by sam on 27/01/2022.
//
/**
 * This is an implementation of a useful fact that the space of linear operators
 * from a vector space U to a vector space V can be identified with the tensor
 * product of U* (the dual of U) tensor product with V. Under this identification
 * the elementary tensor u* tensor v is identified with the operator
 *
 *      x -> <u*, x>v     (x in U).
 *
 * Note that this is well-defined since <u*, x> is an element of the base field
 * so this is simply scalar multiplication of a vector in V. Naturally, most
 * linear operators are represented by arbitrary sums of such elementary
 * tensors. The tensor_defined_operator class implements this kind of operator
 * in terms of a sequence of elementary tensors that are applied sequentially
 * to generate the output.
 *
 */

#ifndef LIBALGEBRA_TENSOR_OPERATOR_H
#define LIBALGEBRA_TENSOR_OPERATOR_H

#include <utility>
#include <vector>

namespace alg {
namespace operators {

template<typename FunctionalType, typename ResultType>
class tensor_defined_operator
{
public:
    using functional_type = FunctionalType;
    using result_type = ResultType;
    using elementary_tensor = std::pair<FunctionalType, ResultType>;

private:
    std::vector<elementary_tensor> m_data;

public:
    explicit tensor_defined_operator(std::initializer_list<elementary_tensor> args) : m_data(args)
    {}

    template <typename Arg>
    result_type operator()(const Arg& arg) const
    {
        result_type result;

        for (const auto& elt : m_data) {
            result.add_scal_prod(elt.second, elt.first(arg));
        }

        return result;
    }
};

}// namespace operators
}// namespace alg
#endif//LIBALGEBRA_TENSOR_OPERATOR_H
