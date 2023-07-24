//
// Created by sam on 08/02/2021.
//

#ifndef LIBALGEBRA_META_H
#define LIBALGEBRA_META_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace alg {
namespace utils {

template<typename T>
struct scalar_base {
private:
    template<typename U, typename R = typename U::SCA>
    static R check(void*);

    template<typename U>
    static U check(...);

public:
    using type = decltype(check<T>(nullptr));
};

template<typename T, typename V>
struct copy_constness {
    typedef V type;
};

template<typename T, typename V>
struct copy_constness<const T, V> {
    typedef typename std::add_const<V>::type type;
};

template<typename T, typename V>
struct copy_constness<const T&, V> {
    typedef typename std::add_const<V>::type type;
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
template<std::size_t... Args>
struct array_holder {
    static constexpr std::array<std::size_t, sizeof...(Args)> data = {Args...};
};

template<std::size_t N, template<std::size_t> class F, std::size_t... Past>
struct generate_array_impl {
    using result = typename generate_array_impl<N - 1, F, F<N>::value, Past...>::result;
};

template<template<std::size_t> class F, std::size_t... Past>
struct generate_array_impl<static_cast<std::size_t>(0), F, Past...> {
    using result = array_holder<F<static_cast<std::size_t>(0)>::value, Past...>;
};

template<std::size_t N, template<std::size_t> class F>
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

template<typename T, template<typename> class MetaFunc, typename>
struct void_or : MetaFunc<T> {
};

template<template<typename> class UNUSED, typename Void>
struct void_or<void, UNUSED, Void>
{
    using type = Void;
};

}// namespace utils
}// namespace alg

#endif// LIBALGEBRA_META_H
