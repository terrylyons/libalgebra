//
// Created by sam on 30/03/2021.
//

#ifndef LIBALGEBRA_ORDER_TRAIT_H
#define LIBALGEBRA_ORDER_TRAIT_H

#include <map>
#include <type_traits>

namespace alg {
namespace utils {

/// Trait indicating whether a map type is ordered
template <typename Map>
struct is_ordered : std::false_type
{};

template <typename K, typename V>
struct is_ordered<std::map<K, V>> : std::true_type
{};

}// namespace utils
}// namespace alg

#endif// LIBALGEBRA_ORDER_TRAIT_H
