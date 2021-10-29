//
// Created by sam on 29/10/2021.
//

#ifndef LIBALGEBRA_OPERATORS_H
#define LIBALGEBRA_OPERATORS_H

namespace alg {
namespace operators {

template<typename Impl, typename ArgumentType, typename ResultType>
class linear_operator : protected Impl
{
public:
    using argument_type = ArgumentType;
    using result_type = ResultType;

protected:
    using implementation_type = Impl;
};

}// namespace operators
}// namespace alg

#endif//LIBALGEBRA_OPERATORS_H
