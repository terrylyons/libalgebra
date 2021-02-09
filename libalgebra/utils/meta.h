//
// Created by sam on 08/02/2021.
//

#ifndef LIBALGEBRA_META_H
#define LIBALGEBRA_META_H

#include <type_traits>

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






} // namespace utils
} // namespace alg


#endif //LIBALGEBRA_META_H
