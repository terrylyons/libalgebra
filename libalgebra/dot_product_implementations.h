//
// Created by sam on 15/11/2021.
//

#ifndef LIBALGEBRA_DOT_PRODUCT_IMPLEMENTATIONS_H
#define LIBALGEBRA_DOT_PRODUCT_IMPLEMENTATIONS_H

#include <libalgebra/vectors.h>

namespace alg {
namespace operators {

template<typename Functional, typename Argument>
class dot_product_implementation
{
public:
    using functional_type = Functional;
    using argument_type = Argument;
    using coefficient_type = typename Functional::coefficient_field;
    using result_type = typename coefficient_type::S;

    template<typename... Args>
    explicit dot_product_implementation(Args&&... args) : m_functional(std::forward<Args>(args)...)
    {}

private:
    using functional_basis_t = typename functional_type::BASIS;
    using argument_basis_t = typename argument_type::BASIS;

    template<template<typename, typename, typename...> class VType, typename... VArgs>
    using fvector_t = vectors::vector<functional_basis_t, coefficient_type, VType, VArgs...>;

    template<template<typename, typename, typename...> class VType, typename... VArgs>
    using avector_t = vectors::vector<argument_basis_t, coefficient_type, VType, VArgs...>;

    template<typename F>
    struct implementation;

    template<template<typename, typename, typename...> class VType, typename... VArgs>
    struct implementation<fvector_t<VType, VArgs...>> {
        template<template<typename, typename, typename...> class AVType, typename... AVArgs>
        static result_type
        eval(const fvector_t<VType, VArgs...>& functional, const avector_t<AVType, AVArgs...>& argument) noexcept
        {
            result_type result(coefficient_type::zero);
            for (auto it : functional) {
                coefficient_type::add_inplace(result, coefficient_type::mul(it.value(), argument[it.key()]));
            }
            return result;
        }
    };

    template<typename... FArgs>
    struct implementation<fvector_t<vectors::dense_vector, FArgs...>> {

        using f_type = fvector_t<vectors::dense_vector, FArgs...>;

        template<template<typename, typename, typename...> class AVType, typename... AVArgs>
        static result_type
        eval(const f_type& functional, const avector_t<AVType, AVArgs...>& argument) noexcept
        {
            result_type result(coefficient_type::zero);
            for (auto it : functional) {
                coefficient_type::add_inplace(result, coefficient_type::mul(it.value(), argument[it.key()]));
            }
            return result;
        }

        template<typename... AArgs>
        static result_type eval(const f_type& functional, const avector_t<vectors::dense_vector, AArgs...>& argument) noexcept
        {
            const auto& dense_functional = functional.base_vector();
            const auto& dense_argument = argument.base_vector();

            auto dim = std::min(dense_functional.dimension(), dense_argument.dimension());

            result_type result(coefficient_type::zero);
            for (DIMN i = 0; i < dim; ++i) {
                coefficient_type::add_inplace(result, coefficient_type::mul(dense_functional.value(i), dense_argument.value(i)));
            }
            return result;
        }
    };

    using impl = implementation<typename Functional::VECT>;

public:
    result_type operator()(const argument_type& argument) const noexcept
    {
        return impl::eval(m_functional, argument);
    }

private:
    functional_type m_functional;
};

}// namespace operators
}// namespace alg

#endif//LIBALGEBRA_DOT_PRODUCT_IMPLEMENTATIONS_H
