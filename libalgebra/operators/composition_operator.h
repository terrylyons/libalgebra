//
// Created by sam on 28/01/2022.
//

#ifndef LIBALGEBRA_COMPOSITION_OPERATOR_H
#define LIBALGEBRA_COMPOSITION_OPERATOR_H

namespace alg {
namespace operators {

template<typename ImplOuter, typename ImplInner, typename ArgumentType, typename ResultType>
class composition_operator : protected ImplInner, protected ImplOuter
{

public:
    explicit composition_operator(const ImplInner& inner, const ImplOuter& outer)
        : ImplInner(inner), ImplOuter(outer)
    {}

    ResultType operator()(const ArgumentType& arg) const
    {
        return ImplOuter::operator()(ImplInner::operator()(arg));
    }
};

}// namespace operators
}// namespace alg

#endif//LIBALGEBRA_COMPOSITION_OPERATOR_H
