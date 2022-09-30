//
// Created by sam on 09/02/2021.
//

#ifndef LIBALGEBRA_BASIS_H
#define LIBALGEBRA_BASIS_H

#include "implementation_types.h"

#include "tags.h"
#include <type_traits>

namespace alg {
namespace basis {

template<typename Basis1, typename Basis2>
struct related_to : std::false_type {
};

template<typename Basis>
struct basis_traits {
    typedef typename Basis::ordering_tag ordering_tag;
    typedef typename Basis::degree_tag degree_tag;
};

}// namespace basis
}// namespace alg

#endif// LIBALGEBRA_BASIS_H
