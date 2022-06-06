//
// Created by sam on 09/02/2021.
//

#ifndef LIBALGEBRA_BASIS_H
#define LIBALGEBRA_BASIS_H

#include "libalgebra/implementation_types.h"

#include "libalgebra/basis/tags.h"
#include <type_traits>

namespace alg {
namespace basis {

template<typename Basis>
struct basis_traits {
    typedef typename Basis::ordering_tag ordering_tag;
    typedef typename Basis::degree_tag degree_tag;
};


template <typename Basis1, typename Basis2>
struct is_related : std::false_type {};


}// namespace basis
}// namespace alg

#endif// LIBALGEBRA_BASIS_H
