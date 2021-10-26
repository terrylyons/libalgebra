//
// Created by sam on 30/03/2021.
//

#ifndef LIBALGEBRA_ORDER_TRAIT_H
#define LIBALGEBRA_ORDER_TRAIT_H

#include <map>

namespace alg {
namespace utils {

/// Trait indicating whether a map type is ordered
template<typename Map>
struct is_ordered {
    static const bool value = false;
};

template<typename K, typename V>
struct is_ordered<std::map<K, V>> {
    static const bool value = true;
};

}// namespace utils
}// namespace alg

#endif// LIBALGEBRA_ORDER_TRAIT_H
