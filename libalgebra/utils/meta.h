//
// Created by sam on 08/02/2021.
//

#ifndef LIBALGEBRA_META_H
#define LIBALGEBRA_META_H

#include <array>
#include <boost/type_traits.hpp>

namespace alg {
namespace utils {

struct true_type {
    static const bool value = true;
};

struct false_type {
    static const bool value = false;
};

// C++98 does not have enable_if, so define our own
template<bool Cond, typename T = void>
struct enable_if {
};

template<typename T>
struct enable_if<true, T> {
    typedef T type;
};

template<typename T, typename V>
struct copy_constness {
    typedef V type;
};

template<typename T, typename V>
struct copy_constness<const T, V> {
    typedef typename boost::add_const<V>::type type;
};

template<typename T, typename V>
struct copy_constness<const T&, V> {
    typedef typename boost::add_const<V>::type type;
};

template<template<unsigned, unsigned> class Compute, unsigned W, unsigned D>
struct populate_array {
    template<typename Array>
    static inline void fill(Array& arr)
    {
        arr[D] = static_cast<typename Array::value_type>(
                std::min(static_cast<unsigned long long>(Compute<W, D>::value), static_cast<unsigned long long>(std::numeric_limits<typename Array::value_type>::max())));
        populate_array<Compute, W, D - 1>::fill(arr);
    }
};

template<template<unsigned, unsigned> class Compute, unsigned W>
struct populate_array<Compute, W, 0> {
    template<typename Array>
    static inline void fill(Array& arr)
    {
        arr[0] = static_cast<typename Array::value_type>(
                std::min(static_cast<unsigned long long>(Compute<W, 0>::value), static_cast<unsigned long long>(std::numeric_limits<typename Array::value_type>::max())));
    }
};

template<bool Cond, typename T1, typename T2>
struct type_selector {
    typedef T1 type;
};

template<typename T1, typename T2>
struct type_selector<true, T1, T2> {
    typedef T2 type;
};

// This pattern is for filling an array at compiletime based on template arguments
// Credit to https://stackoverflow.com/a/2981617/9225581. Slightly modified from
// the original.
template<DIMN... Args>
struct array_holder {
    static constexpr std::array<DIMN, sizeof...(Args)> data = {Args...};
};

template<DIMN N, template<DIMN> class F, DIMN... Past>
struct generate_array_impl {
    using result = typename generate_array_impl<N - 1, F, F<N>::value, Past...>::result;
};

template<template<DIMN> class F, DIMN... Past>
struct generate_array_impl<static_cast<DIMN>(0), F, Past...> {
    using result = array_holder<F<static_cast<DIMN>(0)>::value, Past...>;
};

template<DIMN N, template<DIMN> class F>
struct generate_array {
    using result = typename generate_array_impl<N, F>::result;
};

template<bool Cond, template<typename...> class T1, template<typename...> class T2>
struct template_selector {
    template<typename... A>
    using type = T1<A...>;
};

template<template<typename...> class T1, template<typename...> class T2>
struct template_selector<true, T1, T2> {
    template<typename... A>
    using type = T2<A...>;
};

}// namespace utils
}// namespace alg

#endif// LIBALGEBRA_META_H
