//
// Created by sam on 06/10/2021.
//

#ifndef LIBALGEBRA_DUAL_PAIRING_H
#define LIBALGEBRA_DUAL_PAIRING_H

#include "implimentation_types.h"
#include "vectors/vector.h"
#include "vectors/dense_vector.h"

namespace alg {


/**
 * @brief A dual pairing between two vector spaces
 *
 * A dual pairing is a pair (U, V) together with a bilinear function (u, v) -> \<u, v\>
 * which takes values in the base field. This identifies U with some subspace of the
 * dual space of V, and the function v -> \<u, v\> for each u defines a linear functional
 * on V. This template class declares that such a relationship exists between two
 * vector spaces, defined by the span of two sets of bases, and provides implementations
 * for the bilinear function for the various different implementations of vectors that
 * can exist.
 *
 * As a minimum, this class must implement a public and const operator() which takes
 * a template functional vector and template vector and return a scalar which is the
 * value of \<u, v\> for the given inputs. This can, of course, be specialised for
 * specific types of input vectors to provide optimised implementations; for example,
 * one might wish to give a specialised implementation for dense vectors that uses
 * index in basis order rather than key for its computation.
 *
 * One will not usually use this class directly, and instead use the apply_functional
 * function.
 *
 * @tparam FunctionalBasis Basis of vector space of functionals
 * @tparam VectorBasis Basis of vector space
 * @tparam Coeffs Coefficient field over which to work
 */
template <typename FunctionalBasis, typename VectorBasis, typename Coeffs>
class dual_pairing;


/**
 * @brief Simple implementation of a dual pairing where the action is effectively a dot product
 *
 * This class provides an implementation for a dual pairing where the action is given by the
 * first using the canonical isomorphism of a vector space with dimension n into R^n and then
 * using the dot product. Consequently, the vector spaces of functionals and vectors are assumed
 * to have equal dimension. Moreover, the ordering on the bases must be compatible.
 *
 * This is intended as a utility for quickly implementing a dual pairing without the need to
 * write the dense implementation by hand. To implement, one needs only create a specialisation
 * of the dual_pairing class for the relevant template arguments, and derive the specialisation
 * from this class publically. This will make the action available in the correct way. (See
 * the implementation in tensor.h for an example.)
 *
 * @tparam FunctionalBasis Basis for the vector space of functionals
 * @tparam VectorBasis Basis for the vector space
 * @tparam Coeffs Coefficient field over which to work
 */
template <typename FunctionalBasis, typename VectorBasis, typename Coeffs>
class dot_product_pairing
{
public:

    using functional_basis  = FunctionalBasis;
    using vector_basis      = VectorBasis;
    using scalar_t          = typename Coeffs::S;

    template <template <typename, typename> class Impl>
    using functional_t      = vectors::vector<functional_basis, Coeffs, Impl<functional_basis, Coeffs> >;

    template <template <typename, typename> class Impl>
    using vector_t          = vectors::vector<vector_basis, Coeffs, Impl<vector_basis, Coeffs> >;

    template <
        template<typename, typename> class FunctionalImpl,
        template<typename, typename> class VectorImpl>
    scalar_t operator()(
            const functional_t<FunctionalImpl>& functional,
            const vector_t<VectorImpl>& vector
            ) const
    {
        scalar_t result(Coeffs::zero);
        for (auto& it : functional) {
            Coeffs::add_inplace(result, Coeffs::mul(it.value(), vector[it.key()]));
        }
        return result;
    }

    scalar_t operator()(
        const functional_t<vectors::dense_vector>& functional,
        const vector_t<vectors::dense_vector>& vector
    ) const
    {
        const vectors::dense_vector<functional_basis, Coeffs>& dense_functional
            = vectors::dtl::vector_base_access::convert(functional);
        const vectors::dense_vector<vector_basis, Coeffs>& dense_vector
            = vectors::dtl::vector_base_access::convert(vector);

        DIMN dim = std::min(dense_functional.dimension(), dense_vector.dimension());

        scalar_t result{Coeffs::zero};
        for (DIMN i=0; i < dim; ++i) {
            Coeffs::add_inplace(result, Coeffs::mul(dense_functional.value(i), dense_vector.value(i)));
        }
        return result;
    }


};


/**
 * @brief Apply a functional (via a dual pairing) to a vector from a vector space.
 *
 * This function evaluates, via a dual pairing, the action of an element of a dual
 * space on an element of a vector space. More precisely, if U and V are vector
 * spaces, represented as the types Functional and Vector, respectively, then this
 * function returns the result of evaluating the pairing \<u, v\> for specific
 * elements u of U and v of V.
 *
 * The types Functional and Vector should implement the same interface as
 * vectors::vector and must have the same coefficient_field type. The actual
 * implementation is deferred to the relevant dual_pairing template class
 * which should provide the right implementation for the different types
 * representing elements from the respective vector spaces.
 *
 * @tparam Functional Type of elements of the vector space of functionals
 * @tparam Vector Type of elements of vector space
 * @param functional Element of space of functionals
 * @param vector Element of vector space
 * @return Scalar value of evaluating the functional at the vector.
 */
template <typename Functional, typename Vector>
std::enable_if<std::is_same<typename Functional::coefficient_field, typename Vector::coefficient_field>::value,
               typename Functional::SCALAR>
apply_functional(const Functional& functional, const Vector& vector) {
    using f_basis_t = typename Functional::BASIS;
    using v_basis_t = typename Vector::BASIS;
    using coeff_t   = typename Functional::coefficient_field;

    dual_pairing<f_basis_t, v_basis_t, coeff_t> op;
    return op(functional, vector);
}



}


#endif //LIBALGEBRA_DUAL_PAIRING_H
