//
// Created by sam on 12/02/2021.
//

#ifndef LIBALGEBRAUNITTESTS_HYBRID_VECTOR_H
#define LIBALGEBRAUNITTESTS_HYBRID_VECTOR_H

#include "libalgebra/vectors/base_vector.h"
#include "libalgebra/vectors/sparse_vector.h"
#include "libalgebra/vectors/dense_vector.h"


namespace alg {
namespace vectors {




class basic_resize_manager
{
public:

};



template<
        typename Basis,
        typename Coeffs,
        typename ResizeManager,
        typename DenseStorage,
        typename SparseMap>
class hybrid_vector : dense_vector<Basis, Coeffs, DenseStorage>,
                      sparse_vector<Basis, Coeffs, SparseMap>,
                      base_vector<Basis, Coeffs>
{
    typedef base_vector<Basis, Coeffs> BASE_VEC;
    typedef dense_vector<Basis, Coeffs, DenseStorage> DENSE;
    typedef sparse_vector<Basis, Coeffs, SparseMap> SPARSE;
    typedef ResizeManager MANAGER;

private:

    // The resize manager is responsible for dictating when the
    // vector should resize it's dense part to contain a larger
    // proportion of the vector. This applies only to cases where
    // resizing is optional. This does not apply when resize is
    // mandatory.
    MANAGER m_resize_manager;

public:

    // Type definitions
    typedef Basis BASIS;
    typedef typename BASIS::KEY KEY;

    typedef typename Coeffs::S SCALAR;
    typedef typename Coeffs::Q RATIONAL;

public:

    // Static variables from base_Vec
    using BASE_VEC::basis;
    using BASE_VEC::mone;
    using BASE_VEC::zero;
    using BASE_VEC::one;
    using BASE_VEC::degree_tag;

public:

    // Constructors

    hybrid_vector(void) : DENSE(), SPARSE()
    {}

    explicit hybrid_vector(const KEY& key, const SCALAR s=one)
        : DENSE(), SPARSE(key, s)
    {}

    hybrid_vector(const hybrid_vector& other)
        : DENSE(other), SPARSE(other)
    {}

private:

    hybrid_vector(DENSE dense_vec, SPARSE sparse_vec)
        : DENSE(dense_vec), SPARSE(sparse_vec)
    {}

protected:

    // Resizing the dense part of the vector

    void resize_dense_to_dimension(const DIMN dim)
    {
        if (dim > DENSE::dimension()) {
            DENSE::resize_to_dimension(dim);
            incorporate_sparse();
        } else {
            incorporate_dense(dim);
            DENSE::resize_to_dimension(dim);
        }
    }

    void resize_dense_to_degree(const DEG deg)
    {
        resize_dense_to_dimension(DENSE::start_of_degree(deg));
    }

    void maybe_size()
    {}


private:

    /// Incorporate the sparse elements that should now be dense
    void incorporate_sparse()
    {}

    /// Incorporate the dense elements that should now be sparse
    void incorporate_dense(const DIMN from_index)
    {}

public:

    // Iterator methods




public:

    // Element access

    const SCALAR& operator[](const KEY& key) const
    {
        DIMN idx;
        if ((idx = DENSE::key_to_index(key)) < DENSE::m_dimension) {
            return DENSE::value(idx);
        } else {
            return SPARSE::operator[](key);
        }
    }

    SCALAR& operator[](const KEY& key)
    {
        DIMN idx;
        if ((idx = DENSE::key_to_index(key)) < DENSE::m_dimension) {
            return DENSE::value(idx);
        } else {
            return SPARSE::operator[](key);
        }
    }

public:

    // Comparison operators

    bool operator==(const hybrid_vector& other) const
    {
        // TODO: implement me
        return false;
    }

    bool operator!=(const hybrid_vector& other) const
    {
        return !operator==(other);
    }


public:

    // Arithmetic operators
    // The first few are operators that are simply applied to each
    // part individually.

    hybrid_vector operator-() const
    {
        return hybrid_vector(
                DENSE::operator-(),
                SPARSE::operator-()
            );
    }

    hybrid_vector& operator*=(const SCALAR s)
    {
        DENSE::operator*=(s);
        SPARSE::operator*=(s);
        return *this;
    }

    hybrid_vector& operator/=(const RATIONAL s)
    {
        DENSE::operator/=(s);
        SPARSE::operator/=(s);
        return *this;
    }


    // The next few operators require some matching

    hybrid_vector& operator+=(const hybrid_vector& other)
    {
        DIMN dim = std::max(
                DENSE::m_dimension,
                other.m_dimension
                );
        resize_dense_to_dimension(dim);

        DENSE::operator+=(other);
        SPARSE::operator+=(other);
        incorporate_sparse();
        return *this;
    }

    hybrid_vector& operator-=(const hybrid_vector& other)
    {
        DIMN dim = std::max(
                DENSE::m_dimension,
                other.m_dimension
        );
        resize_dense_to_dimension(dim);

        DENSE::operator-=(other);
        SPARSE::operator-=(other);
        incorporate_sparse();
        return *this;
    }

    __DECLARE_BINARY_OPERATOR(hybrid_vector, *, *=, SCALAR);
    __DECLARE_BINARY_OPERATOR(hybrid_vector, /, /=, RATIONAL);
    __DECLARE_BINARY_OPERATOR(hybrid_vector, +, +=, hybrid_vector);
    __DECLARE_BINARY_OPERATOR(hybrid_vector, -, -=, hybrid_vector);

public:

    // Fused operations

#define DEFINE_FUSED_OP(NAME, ST, OP1, OP2)                                 \
    hybrid_vector& NAME(const KEY& rhs, const ST s)                         \
    {                                                                       \
        DIMN idx;                                                           \
        if ((idx = DENSE::key_to_index(rhs)) < DENSE::dimension()) {        \
            DENSE::value(idx) OP1 (one OP2 s);                              \
        } else {                                                            \
            SPARSE::operator[](rhs) OP1 (one OP2 s);                        \
        }                                                                   \
        return *this;                                                       \
    }                                                                       \
                                                                            \
    hybrid_vector& NAME(const hybrid_vector& rhs, const ST s)               \
    {                                                                       \
        DIMN dim = std::max(DENSE::dimension(), rhs.dimension());           \
        resize_dense_to_dimension(dim);                                     \
        DENSE:: NAME (rhs, s);                                              \
        SPARSE:: NAME (rhs, s);                                             \
        incorporate_sparse();                                               \
        return *this;                                                       \
    }

    DEFINE_FUSED_OP(add_scal_prod, SCALAR, +=, *);
    DEFINE_FUSED_OP(sub_scal_prod, SCALAR, -=, *);
    DEFINE_FUSED_OP(add_scal_div, RATIONAL, +=, /);
    DEFINE_FUSED_OP(sub_scal_div, RATIONAL, -=, /);

#undef DEFINE_FUSED_OP



public:

    // Norms

    SCALAR NormL1() const
    {
        return DENSE::NormL1() + SPARSE::NormL1();
    }

    SCALAR NormL1(const DEG deg) const
    {
        return DENSE::NormL1(deg) + SPARSE::NormL1(deg);
    }

    SCALAR NormLInf() const
    {
        return std::max(DENSE::NormLInf(), SPARSE::NormLInf());
    }

    SCALAR NormLInf(const DEG deg) const
    {
        return std::max(DENSE::NormLInf(deg), SPARSE::NormLInf(deg));
    }

public:

    // Display

    inline friend std::ostream &operator<<(std::ostream &os,
                                           const hybrid_vector &rhs)
    {

    }

public:

    // Transform methods


};








} // namespace vectors
} // namespace alg





#endif //LIBALGEBRAUNITTESTS_HYBRID_VECTOR_H
