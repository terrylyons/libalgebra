//
// Created by sam on 08/02/2021.
//

#ifndef LIBALGEBRA_META_H
#define LIBALGEBRA_META_H

#include <type_traits>
#include <boost/type_traits.hpp>

namespace alg {
namespace utils {

struct true_type
{
    static const bool value = true;
};

struct false_type
{
    static const bool value = false;
};


// C++98 does not have enable_if, so define our own
template <bool Cond, typename T = void>
struct enable_if {};

template <typename T>
struct enable_if<true, T>
{
    typedef T type;
};


template <typename T, typename V>
struct copy_constness
{
    typedef V type;
};

template <typename T, typename V>
struct copy_constness<const T, V>
{
    typedef typename boost::add_const<V>::type type;
};

template <typename T, typename V>
struct copy_constness<const T &, V>
{
    typedef typename boost::add_const<V>::type type;
};

template <template <unsigned, unsigned> class Compute, unsigned W, unsigned D>
struct populate_array
{
    template <typename Array>
    static inline void fill(Array& arr)
    {
        arr[D] = Compute<W, D>::value;
        populate_array<Compute, W, D-1>::fill(arr);
    }

};

template <template <unsigned, unsigned> class Compute, unsigned W>
struct populate_array<Compute, W, 0>
{
    template <typename Array>
    static inline void fill(Array& arr)
    {
        arr[0] = Compute<W, 0>::value;
    }
};





} // namespace utils
} // namespace alg


#endif //LIBALGEBRA_META_H
