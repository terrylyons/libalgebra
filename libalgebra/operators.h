﻿//
// Created by sam on 29/10/2021.
//

#ifndef LIBALGEBRA_OPERATORS_H
#define LIBALGEBRA_OPERATORS_H

#include <type_traits>
#include <utility>

#include "composition_operator.h"
#include "scalar_multiply_operator.h"
#include "sum_operator.h"
#include "vector_bundle.h"

#include <libalgebra/alg_types.h>

namespace alg {
namespace operators {

template<typename Impl, typename ArgumentType, typename ResultType>
class linear_operator : protected Impl
{
    static_assert(
            std::is_same<
                    decltype(std::declval<Impl>()(std::declval<const ArgumentType&>())),
                    ResultType>::value,
            "implementation class must be callable with a const reference to ArgumentType and return ResultType");

public:
    using argument_type = ArgumentType;

protected:
    using implementation_type = Impl;

public:
    using implementation_type::implementation_type;

    using result_type = ResultType;
    using implementation_type::operator();

    vector_bundle<result_type> operator()(const vector_bundle<argument_type>& arg)
    {
        return {implementation_type::operator()(static_cast<const argument_type&>(arg)),
        implementation_type::operator()(arg.fibre())};
    }

protected:
    explicit linear_operator(Impl&& impl) : Impl(std::forward<Impl>(impl))
    {}

public:
    template<typename OtherImpl>
    linear_operator<
            sum_operator<Impl, OtherImpl, ArgumentType, ResultType>,
            ArgumentType,
            ResultType>
    operator+(const linear_operator<OtherImpl, ArgumentType, ResultType>& other) const
    {
        return {*this, other};
    }

    template<typename OtherImpl, typename NewArgumentType>
    linear_operator<
            composition_operator<Impl, OtherImpl, NewArgumentType, ResultType>,
            NewArgumentType,
            ResultType> friend
    operator*(const linear_operator& outer,
              const linear_operator<OtherImpl, NewArgumentType, ArgumentType>& inner)
    {
        using composition_type = composition_operator<Impl, OtherImpl, NewArgumentType, ResultType>;
        using return_type = linear_operator<
                composition_type,
                NewArgumentType,
                ResultType>;

        return return_type(composition_type(static_cast<const OtherImpl&>(inner)), static_cast<const Impl&>(outer));
    }
};

template<typename Impl, typename ArgumentType>
using linear_functional = linear_operator<Impl, ArgumentType, typename ArgumentType::SCALAR>;

namespace dtl {

template<typename Algebra>
class left_multiplication_operator_impl
{
protected:
    using shuffle_algebra_t = Algebra;

public:
    template<typename... Args>
    explicit left_multiplication_operator_impl(Args&&... args)
        : m_lhs(std::forward<Args>(args)...)
    {}

    explicit left_multiplication_operator_impl(shuffle_algebra_t&& alg) : m_lhs(alg)
    {}

    shuffle_algebra_t operator()(const shuffle_algebra_t& arg) const
    {
        return m_lhs * arg;
    }

private:
    shuffle_algebra_t m_lhs;
};

template<typename Algebra>
class right_multiplication_operator_impl
{
protected:
    using shuffle_algebra_t = Algebra;

public:
    template<typename... Args>
    explicit right_multiplication_operator_impl(Args&&... args)
        : m_rhs(std::forward<Args>(args)...)
    {}

    explicit right_multiplication_operator_impl(shuffle_algebra_t&& alg) : m_rhs(alg)
    {}

    shuffle_algebra_t operator()(const shuffle_algebra_t& arg) const
    {
        return arg * m_rhs;
    }

private:
    shuffle_algebra_t m_rhs;
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

namespace dtl {

template<typename ShuffleAlgebra, typename TensorAlgebra>
class adjoint_of_left_multiplication_operator_impl
{
protected:
    using shuffle_algebra_t = ShuffleAlgebra;
    using tensor_algebra_t = TensorAlgebra;

public:
    template<typename... Args>
    explicit adjoint_of_left_multiplication_operator_impl(Args&&... args)
        : m_lhs(std::forward<Args>(args)...)
    {}

    explicit adjoint_of_left_multiplication_operator_impl(tensor_algebra_t&& alg) : m_lhs(alg)
    {}

    shuffle_algebra_t operator()(const shuffle_algebra_t& arg) const
    {
        shuffle_algebra_t result;
        for (auto& pr : m_lhs) {
            result += shift_down(arg, pr.key()) * pr.value();
        }
        return result;

    }

private:
    tensor_algebra_t m_lhs;

    static shuffle_algebra_t shift_down(const shuffle_algebra_t & sh, typename shuffle_algebra_t::KEY word)
    {
        shuffle_algebra_t result(sh), working;
        while (word.size()) {
            auto letter = word.lparent();
            word = word.rparent();
            for (auto pr : result) {
                if (pr.key().size() > 0 && pr.key().lparent() == letter)
                    working[pr.key().rparent()] = result[pr.key()];
            }
            result.swap(working);
            working.clear();
        }
        return result;
    }
};


}// namespace dtl

template<typename ShuffleAlgebra, typename TensorAlgebra>
using adjoint_of_left_multiplication_operator = linear_operator<
        dtl::adjoint_of_left_multiplication_operator_impl<ShuffleAlgebra, TensorAlgebra>,
        ShuffleAlgebra,
        ShuffleAlgebra>;

}// namespace operators
}// namespace alg

#include "free_extension.h"
#include "tensor_operator.h"

#endif//LIBALGEBRA_OPERATORS_H
