//
// Created by sam on 01/02/2021.
//

#ifndef LIBALGEBRA_VECTORS_H
#define LIBALGEBRA_VECTORS_H

#define __DECLARE_BINARY_OPERATOR(T1, NEWOP, OLDOP, T2) \
    T1 operator NEWOP(const T2& rhs) const              \
    {                                                   \
        T1 result(*this);                               \
        return result OLDOP rhs;                        \
    }

#define __DECLARE_UNARY_OPERATOR(NEWT, NEWOP, OLDOP, OLDT) \
    NEWT operator NEWOP(void) const { return OLDT::operator OLDOP(); }

#include "libalgebra/vectors/dense_vector.h"
#include "libalgebra/vectors/hybrid_vector.h"
#include "libalgebra/vectors/sparse_vector.h"
#include "libalgebra/vectors/vector.h"

namespace alg {
namespace vectors {
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

#undef __DECLARE_UNARY_OPERATOR
#undef __DECLARE_BINARY_OPERATOR

#endif// LIBALGEBRA_VECTORS_H
