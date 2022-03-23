//
// Created by sam on 05/02/2021.
//

#ifndef LIBALGEBRA_BASE_VECTOR_H
#define LIBALGEBRA_BASE_VECTOR_H

#include "libalgebra/basis/basis.h"

namespace alg {
namespace vectors {

namespace dtl {

struct access_type_dense {
};
struct access_type_sparse {
};

template<typename vector>
struct data_access_base;

}// namespace dtl

/**
 * @brief base class for vectors
 *
 * provides static members for useful scalar constants, basis instance,
 * and the degree tag, which is used to dispatch to the correct function
 * in various places in the vector hierarchy.
 *
 * all vector types should inherit from this class.
 *
 * @tparam basis basis type
 * @tparam coeff coefficient field
 */
template<typename Basis, typename Coeff>
class base_vector
{
public:
    typedef typename Coeff::S scalar;
    typedef typename Coeff::Q rational;

    typedef typename alg::basis::basis_traits<Basis> basis_traits;

    static const typename basis_traits::degree_tag degree_tag;

    static Basis basis;
    static const scalar one;
    static const scalar mone;
    static const scalar zero;
};

// initialisation of static members of base_vector

/// static initialisation of the sparse_vector basis.
template<typename basis, typename coeff>
basis base_vector<basis, coeff>::basis;

/// static initialisation of the scalar constant +1.
template<typename basis, typename coeff>
const typename coeff::S base_vector<basis, coeff>::one(+1);

/// static initialisation of the scalar constant 0.
template<typename basis, typename coeff>
const typename coeff::S base_vector<basis, coeff>::zero(0);

/// static initialisation of the scalar constant -1.
template<typename basis, typename coeff>
const typename coeff::S base_vector<basis, coeff>::mone(-1);

template<typename basis, typename coeff>
const typename alg::basis::basis_traits<basis>::degree_tag
        base_vector<basis, coeff>::degree_tag;

}// namespace vectors
}// namespace alg

#endif// LIBALGEBRA_BASE_VECTOR_H
