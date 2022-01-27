//
// Created by sam on 27/01/2022.
//

#ifndef LIBALGEBRA_SUM_OPERATOR_H
#define LIBALGEBRA_SUM_OPERATOR_H

namespace alg {
namespace operators {

template<typename Impl1, typename Impl2, typename ArgumentType, typename ResultType>
class sum_operator : protected Impl1, protected Impl2
{

public:
    explicit sum_operator(const Impl1& impl1, const Impl2& impl2)
        : Impl1(impl1), Impl2(impl2)
    {}

    ResultType operator()(const ArgumentType& arg) const
    {
        return Impl1::operator()(arg) + Impl2::operator()(arg);
    }
};

}// namespace operators
}// namespace alg

#endif//LIBALGEBRA_SUM_OPERATOR_H
