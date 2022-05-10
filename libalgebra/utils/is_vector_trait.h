//
// Created by sam on 24/01/2022.
//

#ifndef LIBALGEBRA_IS_VECTOR_TRAIT_H
#define LIBALGEBRA_IS_VECTOR_TRAIT_H

#include <libalgebra/vectors/vector.h>
#include <type_traits>

namespace alg {
namespace utils {

template<typename T>
struct is_vector_type {
private:

    template <typename B, typename C, template <typename, typename, typename...> class Type, typename... Args>
    static std::true_type check(vectors::vector<B, C, Type, Args...>&);

    static std::false_type check(...);

public:

    static constexpr bool value = decltype(check(std::declval<T>()))::value;

};

}// namespace utils
}// namespace alg

#endif//LIBALGEBRA_IS_VECTOR_TRAIT_H
