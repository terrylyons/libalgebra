//
// Created by sam on 05/02/2021.
//

#ifndef LIBALGEBRA_BASE_VECTOR_H
#define LIBALGEBRA_BASE_VECTOR_H

#include "libalgebra/basis/basis.h"
#include <libalgebra/coefficients/coefficients.h>

#include <type_traits>

namespace alg {
namespace vectors {

/**
 * @brief Base class for vectors
 *
 * Provides static members for useful scalar constants, basis instance,
 * and the degree tag, which is used to dispatch to the correct function
 * in various places in the vector hierarchy.
 *
 * All vector types should inherit from this class.
 *
 * @tparam Basis Basis type
 * @tparam Coeff Coefficient field
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

    typedef typename alg::basis::basis_traits<Basis> BASIS_TRAITS;

    static const typename BASIS_TRAITS::degree_tag degree_tag;

    static BASIS basis;

    /*
     * In the future, these will be removed from base_vector,
     * and we will use the constants from the coefficient ring instead.
     */
    static const SCALAR zero;
    static const typename std::enable_if<coefficients_t::is_unital, SCALAR>::type one;
    static const typename std::enable_if<coefficients_t::is_unital, SCALAR>::type mone;
};

// Initialisation of static members of base_vector

/// Static initialisation of the sparse_vector basis.
template<typename Basis, typename Coeff>
Basis base_vector<Basis, Coeff>::basis;

/// Static initialisation of the scalar constant 0.
template<typename Basis, typename Coeff>
const typename Coeff::S base_vector<Basis, Coeff>::zero = Coeff::zero;

/// Static initialisation of the scalar constant +1.
template<typename Basis, typename Coeff>
const typename std::enable_if<Coeff::is_unital, typename Coeff::S>::type base_vector<Basis, Coeff>::one = Coeff::one;

/// Static initialisation of the scalar constant -1.
template<typename Basis, typename Coeff>
const typename std::enable_if<Coeff::is_unital, typename Coeff::S>::type base_vector<Basis, Coeff>::mone = Coeff::mone;

template<typename Basis, typename Coeff>
const typename alg::basis::basis_traits<Basis>::degree_tag
        base_vector<Basis, Coeff>::degree_tag;

}// namespace vectors
}// namespace alg

#endif// LIBALGEBRA_BASE_VECTOR_H
