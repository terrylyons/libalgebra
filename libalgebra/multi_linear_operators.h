//
// Created by sam on 21/01/2022.
//

#ifndef LIBALGEBRA_MULTI_LINEAR_OPERATORS_H
#define LIBALGEBRA_MULTI_LINEAR_OPERATORS_H

#include <libalgebra/vector_bundle.h>

#include <type_traits>
#include <utility>

namespace alg {
namespace operators {

template<typename Impl, typename... ArgumentTypes>
class multi_linear_operator : protected Impl
{
    using result_type = decltype(std::declval<Impl>()(std::declval<const ArgumentTypes&>()...));

    static_assert(
            !std::is_same<
                    result_type,
                    void
            >::value,
            "implementation must be callable with const references to ArgumentTypes and non-void result");

protected:
    using implementation_type = Impl;

public:
    /// Inherit the constructors from the implementation
    using Impl::Impl;

    /// Inherit the operation of the multi-linear map from the implementation
    using implementation_type::operator();


};

}//namespace operators
}//namespace alg

#endif//LIBALGEBRA_MULTI_LINEAR_OPERATORS_H
