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
struct with_degree
{
    static const typename alg::utils::enable_if<(D > 0), DEG>::type
        max_degree = D;
};

struct without_degree {};


template <typename OrderOperator>
struct ordered
{
    typedef OrderOperator order;
};

struct unordered {};



} // namespace basis
} // namespace alg





#endif //LIBALGEBRA_TAGS_H
