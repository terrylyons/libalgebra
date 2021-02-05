//
// Created by sam on 05/02/2021.
//

#ifndef LIBALGEBRA_BASE_VECTOR_H
#define LIBALGEBRA_BASE_VECTOR_H

namespace alg {
namespace vectors {

template<typename Basis, typename Field>
class base_vector {
public:
    typedef typename Basis BASIS;
    typedef typename Field::S SCALAR;
    typedef typename Field::Q RATIONAL;

    static BASIS basis;
    static const SCALAR one;
    static const SCALAR mone;
    static const SCALAR zero;
};


// Initialisation of static members of base_vector

/// Static initialisation of the sparse_vector basis.
template<typename Basis, typename Field>
Basis base_vector<Basis, Field>::basis;

/// Static initialisation of the scalar constant +1.
template<typename Basis, typename Field>
const typename Field::S base_vector<Basis, Field>::one(+1);

/// Static initialisation of the scalar constant 0.
template<typename Basis, typename Field>
const typename Field::S base_vector<Basis, Field>::zero(0);

/// Static initialisation of the scalar constant -1.
template<typename Basis, typename Field>
const typename Field::S base_vector<Basis, Field>::mone(-1);

} // namespace alg
} // namespace vectors

#endif //LIBALGEBRA_BASE_VECTOR_H
