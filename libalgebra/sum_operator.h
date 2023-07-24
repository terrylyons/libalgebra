//
// Created by sam on 27/01/2022.
//

#ifndef LIBALGEBRA_SUM_OPERATOR_H
#define LIBALGEBRA_SUM_OPERATOR_H

namespace alg {
namespace operators {

template<typename Impl1, typename Impl2>
class sum_operator
{
    Impl1 left;
    Impl2 right;

public:
    explicit sum_operator(const Impl1& impl1, const Impl2& impl2)
        : left(impl1), right(impl2)
    {}

    template <typename ArgumentType>
    auto operator()(const ArgumentType& arg) const -> decltype(left(arg) + right(arg))
    {
        return left(arg) + right(arg);
    }
};

}// namespace operators
}// namespace alg

#endif//LIBALGEBRA_SUM_OPERATOR_H
