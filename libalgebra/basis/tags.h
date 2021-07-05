//
// Created by sam on 09/02/2021.
//

#ifndef LIBALGEBRA_TAGS_H
#define LIBALGEBRA_TAGS_H

#include <functional>

#include "libalgebra/utils/meta.h"

namespace alg {
namespace basis {


template <DEG D>
struct with_degree {
    static const typename alg::utils::enable_if<(D > 0),
                                                DEG>::type max_degree;
};

template <DEG D> const typename alg::utils::enable_if<(D > 0),
                                                      DEG>::type with_degree<D>::max_degree = D;


struct without_degree {
};


template <typename OrderOperator>
struct ordered {
    typedef OrderOperator order;

    struct pair_order {
        template <typename Key,
                typename Scalar>
        bool operator()(const std::pair<Key,
                                        Scalar> &p1, const std::pair<Key,
                                                                     Scalar> &p2) const {
            order o;
            return o(p1.first, p2.first);
        }
    };

};

struct unordered {
};


} // namespace basis
} // namespace alg





#endif //LIBALGEBRA_TAGS_H
