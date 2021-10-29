//
// Created by sam on 29/10/2021.
//

#ifndef LIBALGEBRA_OPERATORS_H
#define LIBALGEBRA_OPERATORS_H

#include <type_traits>

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

public:
    using implementation_type::implementation_type;

    using implementation_type::operator();
};

template<typename Impl, typename ArgumentType>
using linear_functional = linear_operator<Impl, ArgumentType, typename ArgumentType::SCALAR>;

namespace dtl {

template<typename Algebra>
class left_multiplication_operator_impl
{
protected:
    typename algebra_t = Algebra;

public:
    template<typename... Args>
    explicit left_multiplication_operator_impl(Args&&... args)
        : m_lhs(std::forward<Args>(args)...)
    {}

    explicit left_multiplication_operator_impl(algebra_t&& alg) : m_lhs(alg)
    {}

    algebra_t operator()(const algebra_t& arg) const
    {
        return m_lhs * arg;
    }

private:
    algebra_t m_lhs;
};

template<typename Algebra>
class right_multiplication_operator_impl
{
protected:
    typename algebra_t = Algebra;

public:
    template<typename... Args>
    explicit right_multiplication_operator_impl(Args&&... args)
        : m_rhs(std::forward<Args>(args)...)
    {}

    explicit right_multiplication_operator_impl(algebra_t&& alg) : m_rhs(alg)
    {}

    algebra_t operator()(const algebra_t& arg) const
    {
        return arg * m_rhs;
    }

private:
    algebra_t m_rhs;
};

}//namespace dtl

template<typename Algebra>
using left_multiplication_operator = linear_operator<
        dtl::left_multiplication_operator_impl<Algebra>,
        Algebra,
        Algebra>;

template<typename Algebra>
using right_multiplication_operator = linear_operator<
        dtl::right_multiplication_operator_impl<Algebra>,
        Algebra,
        Algebra>;

}// namespace operators
}// namespace alg

#endif//LIBALGEBRA_OPERATORS_H
