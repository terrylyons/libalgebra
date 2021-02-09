//
// Created by sam on 08/02/2021.
//
#include <vector>
#include <utility>

#include "libalgebra/vectors/base_vector.h"
#include "libalgebra/utils/meta.h"
#include "libalgebra/basis/basis.h"

#ifndef LIBALGEBRA_DENSE_VECTOR_H
#define LIBALGEBRA_DENSE_VECTOR_H


namespace alg {
namespace vectors {

namespace dtl {

template<typename Basis>
struct requires_order
{
    typedef typename alg::basis::basis_traits<Basis>::ordering_tag::order
            key_ordering;
};

} // namespace dtl

template<typename Basis, typename Coeffs,
        typename Storage = std::vector<typename Coeffs::S>>
class dense_vector : base_vector<Basis, Coeffs>, dtl::requires_order<Basis>
{
    typedef Storage STORAGE;
    typedef base_vector <Basis, Coeffs> BASE_VEC;
    STORAGE m_data;
    DIMN m_dimension;
    DEG m_degree = 0;

public:

    // Type definitions
    typedef Basis BASIS;
    typedef Coeffs COEFFS;
    typedef typename BASIS::KEY KEY;
    typedef typename COEFFS::S SCALAR;
    typedef typename COEFFS::Q RATIONAL;
    typedef typename STORAGE::iterator iterator;
    typedef typename STORAGE::const_iterator const_iterator;

    typedef typename dtl::requires_order<Basis>::key_ordering key_ordering;

public:
    // static variables defined in base_vector
    using BASE_VEC::basis;
    using BASE_VEC::one;
    using BASE_VEC::mone;
    using BASE_VEC::zero;

    using BASE_VEC::degree_tag;

public:

    // Constructors
    /// Default constructor
    dense_vector(void) : m_data(), m_dimension(0)
    {}

    /// Copy constructor
    dense_vector(const dense_vector &v)
            : m_data(v.m_data),
              m_dimension(0)
    {}

    /// Unidimensional constructor
    explicit dense_vector(const KEY &k, const SCALAR &s = one)
            : m_data(), m_dimension(0)
    {
        DIMN idx = resize_for_key(k, degree_tag);
        m_data[idx] = s;
    }

private:

    template<DEG D>
    DIMN resize_for_key(const KEY &key, alg::basis::with_degree<D>)
    {
        DEG d = basis.degree(key);
        resize_to_degree(d);
        return index_from_key(key);
    }

    DIMN resize_for_key(const KEY &key, alg::basis::without_degree)
    {
        DIMN idx = key_to_index(key);
        resize_to_dimension(idx);
        return idx;
    }

    template<DEG D>
    DIMN adjust_dimension(const DIMN dim, alg::basis::with_degree<D>)
    {
        DEG d = index_to_degree(dim);
        return start_of_degree(d + 1);
    }

    DIMN adjust_dimension(const DIMN dim, alg::basis::without_degree)
    {
        return dim;
    }

    template<DEG D>
    void set_degree(alg::basis::with_degree<D>)
    {
        if (m_dimension == 0) {
            m_degree = 0;
        } else {
            m_degree = index_to_degree(m_dimension - 1);
        }
    }

    void set_degree(alg::basis::without_degree)
    {}


public:

    // resizing methods

    /// Reserve to dimension
    void reserve_to_dimension(const DIMN dim)
    {
        m_data.reserve(adjust_dimension(dim, degree_tag));
    }

    /// Reserve to degree
    void
    reserve_to_degree(const DEG deg)
    {
        m_data.reserve(start_of_degree(deg + 1));
    }

    /// Resize to dimension
    void resize_to_dimension(const DIMN dim)
    {
        DIMN new_dim = adjust_dimension(dim, degree_tag);
        m_data.resize(new_dim);
        m_dimension = new_dim;
        set_degree(degree_tag);
    }

    /// Reserve to degree
    void
    resize_to_degree(const DEG deg)
    {
        DEG target_deg = std::min(degree_tag.max_degree, deg);
        DIMN dim = start_of_degree(target_deg + 1);
        m_data.resize(dim);
        m_dimension = dim;
        m_degree = target_deg;
    }


public:

    // iterator methods

    iterator begin()
    { return m_data.begin(); }

    iterator end()
    { return m_data.end(); }

    const_iterator begin() const
    { return m_data.begin(); }

    const_iterator end() const
    { return m_data.end(); }

    const_iterator cbegin() const
    { return m_data.begin(); }

    const_iterator cend() const
    { return m_data.end(); }

    std::pair<iterator, bool> insert(iterator it, SCALAR val)
    {
        //TODO:: implement me
    }

    void insert(const KEY &key, SCALAR val)
    {
        //TODO: implement me
    }

    template<typename InputIterator>
    void insert(InputIterator begin, InputIterator end)
    {
        //TODO: implement me
    }

    void erase(const KEY &key)
    {
        //TODO: implement me
    }

    void erase(const iterator &it)
    {
        //TODO: implement me
    }

public:

    void clear()
    {
        m_data.clear();
        m_dimension = m_data.size();
        set_degree(degree_tag);
    }

public:

    /// Swap operation
    void swap(dense_vector &other)
    {
        m_data.swap(other.m_data);
        std::swap(m_dimension, other.m_dimension);
        std::swap(m_degree, other.m_degree);
    }

public:

    // Access to information about the vector and coefficients

    DIMN dimension() const
    { return m_dimension; }

    DEG degree() const
    { return m_degree; }

    SCALAR &value(const DIMN dim)
    {
        return m_data[dim];
    }

    const SCALAR &value(const DIMN dim) const
    {
        if (dim >= m_dimension) {
            return zero;
        }
        return m_data[dim];
    }

    DIMN size() const
    {
        DIMN sz = 0;
        for (DIMN i = 0; i < m_dimension; ++i) {
            sz += (m_data[i] > 0 ? 1 : 0);
        }
        return sz;
    }

    bool empty() const
    {
        m_data.empty();
    }

private:

    // Index of key and key of index
    static KEY index_to_key(const DIMN index)
    {
        return basis.index_to_key(index);
    }

    static DIMN key_to_index(const KEY &key)
    {
        return basis.key_to_index(key);
    }

    static DEG index_to_degree(const DIMN idx)
    {
        return basis.degree(index_to_key(idx));
    }

    static DIMN start_of_degree(const DEG deg)
    {
        return basis.start_of_degree(deg);
    }

    template <DEG D>
    static DIMN max_dimension(alg::basis::with_degree<D>)
    {
        return basis.start_of_degree(degree_tag.max_degree + 1);
    }

    static DIMN max_dimension(alg::basis::without_degree)
    {
        return 0;
    }

public:

    // Element access (map API)

    const SCALAR &operator[](const KEY &k) const
    {
        return value(key_to_index(k));
    }

    SCALAR &operator[](const KEY &k)
    {
        return value(key_to_index(k));
    }


public:

    // Arithmetic operation

    /// Unary negation
    dense_vector operator-(void) const
    {
        dense_vector result;
        result.resize_to_dimension(m_dimension);

        for (DIMN i = 0; i < m_dimension; ++i) {
            result.m_data[i] = -m_data[i];
        }

        return result;
    }

    /// Inplace scalar multiplication
    dense_vector &operator*=(const SCALAR s)
    {
        if (s == zero) {
            m_data.clear();
        } else {

            for (DIMN i = 0; i < m_dimension; ++i) {
                m_data[i] *= s;
            }
        }
        return *this;
    }

    /// Inplace rational division
    dense_vector &operator/=(const RATIONAL s)
    {
        // Instead of using the /= operator, us the *=
        // The compiler probably applies this optimisation
        // automatically, but let's be safe.
        SCALAR val = one / s;
        for (DIMN i = 0; i < m_dimension; ++i) {
            m_data[i] *= val;
        }
        return *this;
    }

    /// Scalar multiplication
    __DECLARE_BINARY_OPERATOR(dense_vector, *
    , *=, SCALAR);

    /// Rational division
    __DECLARE_BINARY_OPERATOR(dense_vector,
    /, /=, RATIONAL);

    /// Inplace addition
    dense_vector &operator+=(const dense_vector &rhs)
    {
        if (empty())
            return (*this = rhs);

        if (rhs.empty())
            return *this;

        if (rhs.m_dimension > m_dimension) {
            resize_to_dimension(rhs.m_dimension);
        }

        for (DIMN i = 0; i < rhs.m_dimension; ++i) {
            m_data[i] += rhs.m_data[i];
        }

        return *this;
    }

    /// Inplace subtraction
    dense_vector &operator-=(const dense_vector &rhs)
    {
        if (rhs.empty())
            return *this;

        if (rhs.m_dimension > m_dimension) {
            resize_to_dimension(rhs.m_dimension);
        }

        for (DIMN i = 0; i < rhs.m_dimension; ++i) {
            m_data[i] -= rhs.m_data[i];
        }

        return *this;
    }

    /// Addition
    __DECLARE_BINARY_OPERATOR(dense_vector,
    +, +=, dense_vector);

    /// Subtraction
    __DECLARE_BINARY_OPERATOR(dense_vector,
    -, -=, dense_vector);

    /// Inplace coordinatewise min
    dense_vector &operator&=(const dense_vector &rhs)
    {
        DIMN mid = std::min(m_dimension, rhs.m_dimension);

        for (DIMN i = 0; i < mid; ++i) {
            m_data[i] = std::min(m_data[i], rhs.m_data[i]);
        }

        if (m_dimension < rhs.m_dimension) {
            resize_to_dimension(rhs.m_dimension);

            for (DIMN i = mid; i < rhs.m_dimension; ++i) {
                m_data[i] = std::min(zero, rhs.m_data[i]);
            }
        } else if (m_dimension > rhs.m_dimension) {

            for (DIMN i = mid; i < m_dimension; ++i) {
                m_data[i] = std::min(zero, rhs.m_data[i]);
            }
        }
        return *this;
    }

    /// Inplace coordinatewise max
    dense_vector &operator|=(const dense_vector &rhs)
    {
        DIMN mid = std::min(m_dimension, rhs.m_dimension);

        for (DIMN i = 0; i < mid; ++i) {
            m_data[i] = std::max(m_data[i], rhs.m_data[i]);
        }

        if (m_dimension < rhs.m_dimension) {
            resize_to_dimension(rhs.m_dimension);

            for (DIMN i = mid; i < rhs.m_dimension; ++i) {
                m_data[i] = std::max(zero, rhs.m_data[i]);
            }
        }

        for (DIMN i = mid; i < m_dimension; ++i) {
            m_data[i] = std::max(zero, m_data[i]);
        }

        return *this;
    }

    /// Coordinatewise min
    __DECLARE_BINARY_OPERATOR(dense_vector, &
    , &=, dense_vector);

    /// Coordinatewise max
    __DECLARE_BINARY_OPERATOR(dense_vector,
    |, |=, dense_vector);


public:

#define DECLARE_FUSED_OP(NAME, OP1, OP2, T)                           \
    dense_vector& NAME(const KEY& rhs, const T s) {                   \
        operator[](rhs) OP1 (one OP2 s);                              \
        return *this;                                                 \
    }                                                                 \
                                                                      \
    dense_vector& NAME(const dense_vector& rhs, const T s) {          \
        if (m_dimension < rhs.m_dimension) {                          \
            resize_to_dimension(rhs.m_dimension);                     \
        }                                                             \
                                                                      \
        for (DIMN i=0; i<rhs.m_dimension; ++i) {                      \
            m_data[i] OP1 (rhs.m_data[i] OP2 s);                      \
        }                                                             \
        return *this;                                                 \
    }

    // Fused add/sub scalar multiplication/division
    DECLARE_FUSED_OP(add_scal_prod, +=, *, SCALAR);

    DECLARE_FUSED_OP(sub_scal_prod, -=, *, SCALAR);

    DECLARE_FUSED_OP(add_scal_div, +=, *, RATIONAL);

    DECLARE_FUSED_OP(sub_scal_div, +=, *, RATIONAL);

#undef DECLARE_FUSED_OP

public:

    bool operator==(const dense_vector &rhs) const
    {
        DIMN mid = std::min(m_dimension, rhs.m_dimension);

        for (DIMN i = 0; i < mid; ++i) {
            if (m_data[i] != rhs.m_data[i]) {
                return false;
            }
        }

        for (DIMN i = mid; i < m_dimension; ++i) {
            if (m_data[i] != zero) {
                return false;
            }
        }

        for (DIMN i = mid; i < rhs.m_dimension; ++i) {
            if (rhs.m_data[i] != zero) {
                return false;
            }
        }
        return true;
    }

    bool operator!=(const dense_vector &rhs) const
    {
        return !operator==(rhs);
    }

public:

    // Norms

    SCALAR NormL1() const
    {
        SCALAR ans(zero);
        for (DIMN i = 0; i < m_dimension; ++i) {
            ans += abs(m_data[i]);
        }
        return ans;
    }

    SCALAR NormL1(const DEG deg) const
    {
        SCALAR ans(zero);
        for (DIMN i = start_of_degree(deg); i < start_of_degree(deg + 1); ++i) {
            ans += abs(m_data[i]);
        }
        return ans;
    }

    SCALAR NormLInf() const
    {
        SCALAR ans(zero);
        for (DIMN i = 0; i < m_dimension; ++i) {
            SCALAR abs_val = abs(m_data[i]);
            ans = (abs_val > ans) ? abs_val : ans;
        }
    }

    SCALAR NormLInf(const DEG deg) const
    {
        SCALAR ans(zero);
        for (DIMN i = start_of_degree(deg); i < start_of_degree(deg + 1); ++i) {
            SCALAR abs_val = abs(m_data[i]);
            ans = (abs_val > ans) ? abs_val : ans;
        }
        return ans;
    }

public:

    static bool comp(typename std::pair<KEY, SCALAR> lhs,
                     typename std::pair<KEY, SCALAR> rhs)
    {
        return lhs.first < rhs.first;
    }

public:

    inline friend std::ostream &operator<<(std::ostream &os,
                                           const dense_vector &rhs)
    {
        std::pair < BASIS * , KEY > token;
        token.first = &dense_vector::basis;

        os << '{';
        for (DIMN i = 0; i < rhs.m_dimension; ++i) {
            token.second = index_to_key(i);
            os << ' ' << rhs.m_data[i] << '(' << token << ')';
        }
        os << ' ' << '}';
        return os;
    }


};

} // namespace vectors
} // namespace alg

#endif //LIBALGEBRA_DENSE_VECTOR_H
