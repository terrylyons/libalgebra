//
// Created by sam on 01/02/2021.
//

#ifndef LIBALGEBRA_VECTORS_H
#define LIBALGEBRA_VECTORS_H

#include "basis.h"
#include "dense_vector.h"
#include "hybrid_vector.h"
#include "sparse_vector.h"
#include "vector.h"

namespace alg {
namespace vectors {

template <typename Vector>
struct vector_traits
{
    using basis_type = typename Vector::BASIS;
    using basis_traits = basis::basis_traits<basis_type>;

    template <typename Basis>
    using is_compatible = basis::related_to<basis_type, Basis>;
};

template<typename Basis, typename Coeffs>
struct vector_type_selector {
    typedef sparse_vector<Basis, Coeffs> sparse_vec;
    typedef dense_vector<Basis, Coeffs> dense_vect;
    typedef hybrid_vector<Basis, Coeffs> hybrid_vect;

    typedef typename alg::utils::type_selector<boost::is_pod<typename Coeffs::S>::value, sparse_vec,
                                               dense_vect>::type type;
};


template <typename Basis, typename Coeffs>
struct template_vector_type_selector
{
    template<typename B, typename C>
    using type = typename alg::utils::template_selector<
            boost::is_pod<typename Coeffs::S>::value,
            sparse_vector,
            dense_vector
    >::template type<B, C>;
};


}// namespace vectors
}// namespace alg


#endif// LIBALGEBRA_VECTORS_H
