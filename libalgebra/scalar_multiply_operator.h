//
// Created by sam on 27/01/2022.
//

#ifndef LIBALGEBRA_SCALAR_MULTIPLY_OPERATOR_H
#define LIBALGEBRA_SCALAR_MULTIPLY_OPERATOR_H

#include <utility>

namespace alg {
namespace operators {

template<typename Impl, typename Scalar>
class scalar_multiply_operator : protected Impl
{
    Scalar m_scalar;

public:
    explicit scalar_multiply_operator(const Impl& implementation, const Scalar& s)
        : Impl(implementation), m_scalar(s)
    {}

    template <typename ArgumentType>
    auto operator()(const ArgumentType& arg) const -> decltype(Impl::operator()(arg)*m_scalar)
    {
        return Impl::operator()(arg) * m_scalar;
    }

    scalar_multiply_operator& operator*=(const Scalar& other)
    {
        m_scalar *= other;
        return *this;
    }

    scalar_multiply_operator operator*(const Scalar& other) const
    {
        return {*this, m_scalar * other};
    }
};

}// namespace operators
}// namespace alg

#endif//LIBALGEBRA_SCALAR_MULTIPLY_OPERATOR_H
