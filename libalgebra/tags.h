//
// Created by sam on 09/02/2021.
//

#ifndef LIBALGEBRA_TAGS_H
#define LIBALGEBRA_TAGS_H

#include <type_traits>
#include <utility>

namespace alg {
namespace basis {

/**
 * @brief Tag to indicate that a basis has an associated degree
 *
 * This tag indicates that the basis has a degree associated with it
 * and communicates the maximum degree values that can be stored.
 *
 * @tparam D Maximum degree, must be non-zero
 */
template<DEG D>
struct with_degree {
    static const typename std::enable_if<(D > 0), DEG>::type max_degree;
};

template<DEG D>
const typename std::enable_if<(D > 0), DEG>::type
        with_degree<D>::max_degree = D;

/**
 * @brief Tag to indicate that a basis has no associated degree
 *
 */
struct without_degree {
};

/**
 * @brief Tag to indicate that the basis as an ordering
 *
 * Indicates that a basis is ordered and provides the ordering for both basis keys
 * and key-value pairs.
 *
 * @tparam OrderOperator Order operator for basis keys
 */
template<typename OrderOperator>
struct ordered {
    typedef OrderOperator order;

    struct pair_order {
        template<typename Key, typename Scalar>
        bool operator()(const std::pair<Key, Scalar>& p1, const std::pair<Key, Scalar>& p2) const
        {
            order o;
            return o(p1.first, p2.first);
        }
    };
};

/// Tag to indicate that a basis has no ordering
struct unordered {
};

}// namespace basis
}// namespace alg

#endif// LIBALGEBRA_TAGS_H
