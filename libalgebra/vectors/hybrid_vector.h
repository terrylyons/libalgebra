//
// Created by sam on 12/02/2021.
//

#ifndef LIBALGEBRA_HYBRID_VECTOR_H
#define LIBALGEBRA_HYBRID_VECTOR_H

#include "libalgebra/vectors/base_vector.h"
#include "libalgebra/vectors/sparse_vector.h"
#include "libalgebra/vectors/dense_vector.h"


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
    typedef base_vector <Basis, Coeffs> BASE_VEC;
    typedef dense_vector <Basis, Coeffs, DenseStorage> DENSE;
    typedef sparse_vector <Basis, Coeffs, SparseMap> SPARSE;
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

    typedef typename dtl::requires_order<Basis>::key_ordering key_ordering;

    typedef void iterator;
    typedef void const_iterator;

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

    explicit hybrid_vector(const KEY &key, const SCALAR s = one)
            : DENSE(), SPARSE(key, s)
    {}

    hybrid_vector(const hybrid_vector &other)
            : DENSE(other), SPARSE(other)
    {}

private:

    hybrid_vector(DENSE dense_vec, SPARSE sparse_vec)
            : DENSE(dense_vec), SPARSE(sparse_vec)
    {}

protected:

    // Resizing the dense part of the vector

    void reserve_dense_to_dimension(const DIMN dim)
    {
        DENSE::resize_to_dimension(dim);
    }

    void reserve_dense_to_degree(const DEG deg)
    {
        DENSE::resize_to_degree(deg);
    }

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
    {
        if (dense_dimension() == 0 || sparse_empty()) {
            return;
        }

        typename SPARSE::iterator it(SPARSE::begin()), end(SPARSE::end());


        DIMN idx;

        while (it != end) {
            if ((idx = key_to_index(it->first)) < dense_dimension()) {
                DENSE::value(idx) += it->second;
                SPARSE::erase(it++);
            } else {
                ++it;
            }
        }

    }

    /// Incorporate the dense elements that will need to be sparse
    void incorporate_dense(const DIMN from_index)
    {}

public:

    // Vector information methods

    DIMN dense_dimension() const
    { return DENSE::dimension(); }

    DIMN max_dense_dimension() const
    {
        return DENSE::max_dimension();
    }

    DEG dense_degree() const
    { return DENSE::degree(); }

    DIMN dense_size() const
    { return DENSE::size(); }

    DIMN sparse_size() const
    { return SPARSE::size(); }

    bool empty() const
    {
        return DENSE::empty() && SPARSE::empty();
    }

    bool dense_empty() const
    { return DENSE::empty(); }

    bool sparse_empty() const
    { return SPARSE::empty(); }

protected:

    // Sparse part and dense part access

    DENSE &dense_part()
    { return *this; }

    const DENSE &dense_part() const
    { return *this; }

    SPARSE &sparse_part()
    { return *this; }

    const SPARSE &sparse_part() const
    { return *this; }

protected:

    // Information about keys from dense vector

    using DENSE::index_to_key;
    using DENSE::key_to_index;

public:

    // Iterator methods

    iterator begin()
    {

    }

    iterator end()
    {}

    const_iterator begin() const
    {}

    const_iterator end() const
    {}

    const_iterator cbegin() const
    {}

    const_iterator cend() const
    {}

    void insert(iterator it, SCALAR val)
    {}

    void insert(const KEY &key, SCALAR val)
    {
        operator[](key) = val;
    }

    template<typename InputIterator>
    void insert(InputIterator begin, InputIterator end)
    {

    }

    const_iterator find(const KEY &key) const
    {}

    iterator find(const KEY &key)
    {}

    void erase(const KEY &key)
    {}

    void erase(iterator &it)
    {}

public:

    void clear()
    {
        DENSE::clear();
        SPARSE::clear();
    }

public:

    void swap(hybrid_vector& other)
    {
        DENSE::swap(other);
        SPARSE::swap(other);
    }

public:

    // Element access

    const SCALAR &operator[](const KEY &key) const
    {
        DIMN idx;
        if ((idx = key_to_index(key)) < DENSE::dimension()) {
            return DENSE::value(idx);
        } else {
            return SPARSE::operator[](key);
        }
    }

    SCALAR &operator[](const KEY &key)
    {
        DIMN idx;
        if ((idx = key_to_index(key)) < DENSE::dimension()) {
            return DENSE::value(idx);
        } else {
            return SPARSE::operator[](key);
        }
    }

    SCALAR& value(const DIMN idx)
    {
        if (idx < DENSE::dimension()) {
            return DENSE::value(idx);
        } else {
            return SPARSE::operator[](index_to_key(idx));
        }
    }

    const SCALAR& value(const DIMN idx) const
    {
        if (idx < DENSE::dimension()) {
            return DENSE::value(idx);
        } else {
            return SPARSE::operator[](index_to_key(idx));
        }
    }

public:

    // Comparison operators

    bool operator==(const hybrid_vector &rhs) const
    {
        if (empty() && rhs.empty()) {
            return true;
        }

        if (dense_dimension() == 0 && rhs.dense_dimension() == 0) {
            return sparse_part() == rhs.sparse_part();
        }

        if (sparse_empty() && rhs.sparse_empty()) {
            return dense_part() == rhs.dense_part();
        }


        std::pair<DIMN, bool> dense_part_eq = DENSE::equal_to_min(rhs);
        if (dense_part_eq.first != 0 && !dense_part_eq.second) {
            return false;
        }

        typename SPARSE::const_iterator cit, cend, oit;

        if (dense_part_eq.first == dense_dimension()
            && dense_part_eq.first == rhs.dense_dimension()) {
            return sparse_part() == rhs.sparse_part();
        } else if (dense_part_eq.first == dense_dimension()) {
            cend = SPARSE::end();

            for (DIMN i = dense_part_eq.first; i < rhs.dense_dimension(); ++i) {
                if (rhs.dense_part().value(i) == zero) {
                    continue;
                }

                cit = SPARSE::find(index_to_key(i));
                if (cit == cend || cit->second != rhs.dense_part().value(i)) {
                    return false;
                }
            }

            for (oit= rhs.sparse_part().begin(); oit != rhs.sparse_part().end(); ++oit) {
                cit = SPARSE::find(oit.first);
                if (cit == cend || cit->second != oit->second) {
                    return false;
                }
            }

            return true;
        } else if (dense_part_eq.first == rhs.dense_dimension()) {

            cend = rhs.sparse_part().end();

            for (DIMN i=dense_part_eq.first; i<dense_dimension(); ++i) {
                if (DENSE::value(i) == zero) {
                    continue;
                }

                cit = rhs.sparse_part().find(index_to_key(i));
                if (cit == cend || cit->second != DENSE::value(i)) {
                    return false;
                }
            }

            for (oit = sparse_part().begin(); oit != sparse_part().end(); ++oit) {
                cit = rhs.sparse_part().find(oit.first);
                if (cit == cend || cit->second != oit->second) {
                    return false;
                }
            }
            return true;

        }

        // Should be unreachable
        assert(false);
        return false;


    }

    bool operator!=(const hybrid_vector &other) const
    {
        return !operator==(other);
    }


public:

    // Arithmetic operators
    // The first few are operators that are simply applied to each
    // part individually.

    hybrid_vector operator-(void) const
    {
        return hybrid_vector(
                DENSE::operator-(),
                SPARSE::operator-()
        );
    }

    hybrid_vector &operator*=(const SCALAR s)
    {
        DENSE::operator*=(s);
        SPARSE::operator*=(s);
        return *this;
    }

    hybrid_vector &operator/=(const RATIONAL s)
    {
        DENSE::operator/=(s);
        SPARSE::operator/=(s);
        return *this;
    }


    // The next few operators require some matching

    hybrid_vector &operator+=(const hybrid_vector &other)
    {
        DIMN dim = std::max(
                DENSE::dimension(),
                other.dimension()
        );
        resize_dense_to_dimension(dim);

        DENSE::operator+=(other);
        SPARSE::operator+=(other);
        incorporate_sparse();
        return *this;
    }

    hybrid_vector &operator-=(const hybrid_vector &other)
    {
        DIMN dim = std::max(
                DENSE::dimension(),
                other.dimension()
        );
        resize_dense_to_dimension(dim);

        DENSE::operator-=(other);
        SPARSE::operator-=(other);
        incorporate_sparse();
        return *this;
    }

    hybrid_vector& operator&=(const hybrid_vector& rhs)
    {
        DIMN dim = std::max(
                DENSE::dimension(),
                rhs.dimension()
        );
        resize_dense_to_dimension(dim);

        DENSE::operator&=(rhs);
        SPARSE::operator&=(rhs);
        incorporate_sparse();
        return *this;
    }

    hybrid_vector& operator|=(const hybrid_vector& rhs)
    {
        DIMN dim = std::max(
                DENSE::dimension(),
                rhs.dimension()
        );
        resize_dense_to_dimension(dim);

        DENSE::operator|=(rhs);
        SPARSE::operator|=(rhs);
        incorporate_sparse();
        return *this;
    }


    __DECLARE_BINARY_OPERATOR(hybrid_vector, *, *=, SCALAR);

    __DECLARE_BINARY_OPERATOR(hybrid_vector, /, /=, RATIONAL);

    __DECLARE_BINARY_OPERATOR(hybrid_vector, +, +=, hybrid_vector);

    __DECLARE_BINARY_OPERATOR(hybrid_vector, -, -=, hybrid_vector);

    __DECLARE_BINARY_OPERATOR(hybrid_vector, &, &=, hybrid_vector);

    __DECLARE_BINARY_OPERATOR(hybrid_vector, |, |=, hybrid_vector);


public:

    // Fused operations

    DEFINE_FUSED_OP(add_scal_prod, SCALAR, +=, *);

    DEFINE_FUSED_OP(sub_scal_prod, SCALAR, -=, *);

    DEFINE_FUSED_OP(add_scal_div, RATIONAL, +=, /);

    DEFINE_FUSED_OP(sub_scal_div, RATIONAL, -=, /);


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
        std::pair<BASIS*, KEY> token;
        token.first = &basis;

        os << '{';

        for (const_iterator cit(rhs.begin()); cit != rhs.end(); ++cit) {
            if (zero != cit->second) {
                token.second = cit->second;
                os << ' ' << cit->second << '(' << token << ')';
            }
        }

        os << ' ' << '}';
        return os;
    }

public:

    static bool comp(std::pair<KEY, SCALAR> lhs,
                     std::pair<KEY, SCALAR> rhs)
    {
        key_ordering ord;
        return ord(lhs.first, rhs.first);
    }


public:

    // Transform methods


};

#undef DEFINE_FUSED_OP

} // namespace vectors
} // namespace alg





#endif //LIBALGEBRA_HYBRID_VECTOR_H
