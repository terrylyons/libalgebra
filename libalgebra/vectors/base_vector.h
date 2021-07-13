//
// Created by sam on 05/02/2021.
//

#ifndef LIBALGEBRA_BASE_VECTOR_H
#define LIBALGEBRA_BASE_VECTOR_H

#include "libalgebra/basis/basis.h"

namespace alg {
namespace vectors {

template <typename Basis, typename Coeff>
class base_vector
{
public:
    typedef Basis BASIS;
    typedef typename Coeff::S SCALAR;
    typedef typename Coeff::Q RATIONAL;

    typedef typename alg::basis::basis_traits<Basis> BASIS_TRAITS;

    static const typename BASIS_TRAITS::degree_tag degree_tag;

    static BASIS basis;
    static const SCALAR one;
    static const SCALAR mone;
    static const SCALAR zero;
};

// Initialisation of static members of base_vector

/// Static initialisation of the sparse_vector basis.
template <typename Basis, typename Coeff> Basis base_vector<Basis, Coeff>::basis;

/// Static initialisation of the scalar constant +1.
template <typename Basis, typename Coeff> const typename Coeff::S base_vector<Basis, Coeff>::one(1);

/// Static initialisation of the scalar constant 0.
template <typename Basis, typename Coeff> const typename Coeff::S base_vector<Basis, Coeff>::zero(0);

/// Static initialisation of the scalar constant -1.
template <typename Basis, typename Coeff> const typename Coeff::S base_vector<Basis, Coeff>::mone(-1);

template <typename Basis, typename Coeff> const typename alg::basis::basis_traits<Basis>::degree_tag
        base_vector<Basis, Coeff>::degree_tag;

} // namespace vectors
} // namespace alg

#endif // LIBALGEBRA_BASE_VECTOR_H
