//
// Created by sam on 28/01/2022.
//

#ifndef LIBALGEBRA_COMPOSITION_OPERATOR_H
#define LIBALGEBRA_COMPOSITION_OPERATOR_H

namespace alg {
namespace operators {

template<typename ImplOuter, typename ImplInner>
class composition_operator : protected ImplInner, protected ImplOuter
{

public:
    explicit composition_operator(const ImplInner& inner, const ImplOuter& outer)
        : ImplInner(inner), ImplOuter(outer)
    {}

    template <typename ArgumentType>
    auto operator()(const ArgumentType& arg) const -> decltype(std::declval<const ImplOuter&>()(std::declval<const ImplInner&>()(arg)))
    {
        return ImplOuter::operator()(ImplInner::operator()(arg));
    }
};

}// namespace operators
}// namespace alg

#endif//LIBALGEBRA_COMPOSITION_OPERATOR_H
