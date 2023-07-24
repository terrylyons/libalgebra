//
// Created by sam on 05/02/2021.
//

#ifndef LIBALGEBRA_BASE_VECTOR_H
#define LIBALGEBRA_BASE_VECTOR_H

#include "libalgebra/basis.h"
#include "libalgebra/coefficients.h"

#include <type_traits>

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
    static_assert(
            coefficients::is_coefficient_ring<Coeff>::value,
            "Coeff should be a coefficient_ring subclass"
            );

public:

    using coefficients_t = Coeff;

    typedef Basis BASIS;
    typedef typename Coeff::S SCALAR;
    typedef typename Coeff::Q RATIONAL;


    typedef typename alg::basis::basis_traits<Basis> basis_traits;

    static const typename basis_traits::degree_tag degree_tag;


    static BASIS basis;

    /*
     * In the future, these will be removed from base_vector,
     * and we will use the constants from the coefficient ring instead.
     */
    static const SCALAR zero;
    static const typename std::enable_if<coefficients_t::is_unital, SCALAR>::type one;
    static const typename std::enable_if<coefficients_t::is_unital, SCALAR>::type mone;

};

// initialisation of static members of base_vector

/// static initialisation of the sparse_vector basis.
template<typename basis, typename coeff>
basis base_vector<basis, coeff>::basis;

/// Static initialisation of the scalar constant 0.
template<typename Basis, typename Coeff>
const typename Coeff::S base_vector<Basis, Coeff>::zero {};

/// Static initialisation of the scalar constant +1.
template<typename Basis, typename Coeff>
const typename std::enable_if<Coeff::is_unital, typename Coeff::S>::type base_vector<Basis, Coeff>::one(1);

/// Static initialisation of the scalar constant -1.
template<typename Basis, typename Coeff>
const typename std::enable_if<Coeff::is_unital, typename Coeff::S>::type base_vector<Basis, Coeff>::mone(-1);


template<typename basis, typename coeff>
const typename alg::basis::basis_traits<basis>::degree_tag
        base_vector<basis, coeff>::degree_tag;

}// namespace vectors
}// namespace alg

#endif// LIBALGEBRA_BASE_VECTOR_H
