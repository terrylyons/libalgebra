//
// Created by sam on 08/02/2021.
//
#include <vector>
#include <utility>

#include "libalgebra/vectors/base_vector.h"
#include "libalgebra/utils/meta.h"

#ifndef LIBALGEBRA_DENSE_VECTOR_H
#define LIBALGEBRA_DENSE_VECTOR_H

template <typename Basis>
struct has_degree
{
    static const bool value = Basis::MAX_DEGREE > 0;
};



typedef std::size_t DIMN;


namespace alg {
namespace vectors {

template<typename Basis, typename Coeffs,
        typename Storage = std::vector<typename Coeffs::S>>
class dense_vector : base_vector<Basis, Coeffs> {

    typedef Storage STORAGE;
    typedef base_vector<Basis, Coeffs> BASE_VEC;

    STORAGE m_data;
    DIMN m_dimension;
    DEG m_degree = 0;

public:

    // Type defintions
    typedef Basis BASIS;
    typedef Coeffs COEFFS;
    typedef typename BASIS::KEY KEY;
    typedef typename COEFFS:S SCALAR;
    typedef typename COEFFS::Q RATIONAL;
    typedef typename STORAGE::iterator iterator;
    typedef typename STORAGE::const_iterator const_iterator;

public:
    // static variables defined in base_vector
    using BASE_VEC::basis;
    using BASE_VEC::one;
    using BASE_VEC::mone;
    using BASE_VEC::zero;

public:

    // Constructors
    /// Default constructor
    dense_vector(void) : STORAGE() {}

    /// Copy constructor
    dense_vector(const dense_vector& v) : STORAGE(v.m_data) {}

    /// Unidimensional constructor
    explicit dense_vector(const KEY& k, const SCALAR& = one) : STORAGE()
    {
        //TODO: Implement me
    }

public:

    // resizing methods

    /// Reserve to dimension
    void reserve_to_dimension(const DIMN dim)
    {
    }

    /// Reserve to degree
    typename alg::utils::enable_if<has_degree<BASIS>::value>::type
    reserve_to_degree(const DEG deg)
    {}

    /// Resize to dimension
    void resize_to_dimension(const DIMN dim)
    {
    }

    /// Reserve to degree
    typename alg::utils::enable_if<has_degree<BASIS>::value>::type
    reserve_to_degree(const DEG deg)
    {}


public:

    // iterator methods

    iterator begin() { return m_data.begin(); }
    iterator end() { return m_data.end(); }
    const_iterator begin() const { return m_data.begin(); }
    const_iterator end() const { return m_data.end(); }

    std::pair<iterator, bool> insert(iterator it, SCALAR val)
    {}

    void insert(const KEY& key, SCALAR val)
    {
    }

    template <typename InputIterator>
    void insert(InputIterator begin, InputIterator end)
    {}

    void erase(const KEY& key)
    {}

    void erase(const iterator& it)
    {}

public:

    void clear()
    {
        m_data.clear();
        m_dimension = m_data.size();
        m_degree = 0;
    }

public:

    /// Swap operation
    void swap(dense_vector& other)
    {
        m_data.swap(other.m_data);
        m_dimension.swap(other.m_dimension);
        m_degree.swap(other.m_degree);
    }

public:

    // Access to information about the vector and coefficients

    DIMN dimension() const { return m_dimension; }

    DEG degree() const { return m_degree; }

    SCALAR& value(const DIMN dim) { return m_data[dim]; }
    const SCALAR& value(const DIMN dim) const { return m_data[dim]; }

    DIMN size() const
    {
        DIMN sz = 0;
        for (DIMN i=0; i<m_dimension; ++i) {
            sz += (m_data[i] > 0 ? 1 : 0);
        }
        return sz;
    }

    bool empty() const
    {
        return m_data.empty();
    }

private:

    // Index of key and key of index
    static KEY key_from_index(const DIMN index)
    {
        //TODO: implement me
        return 0;
    }

    static DIMN index_from_key(const KEY& key)
    {
        //TODO: implement me
        return KEY();
    }



public:

    // Element access (map API)

    const SCALAR& operator[](const KEY& k) const
    {
        return value(index_from_key(k));
    }

    SCALAR& operator[](const KEY& k)
    {
        return value(index_from_key(k));
    }


public:

    // Arithmetic operation

    /// Unary negation
    dense_vector operator-(void) const
    {
        dense_vector result;
        result.resize_to_dimension(m_dimension);

        for (DIMN i=0; i<m_dimension; ++i) {
            result.m_data[i] = -m_data[i];
        }

        return result;
    }

    /// Inplace scalar multiplication
    dense_vector& operator*=(const SCALAR s)
    {
        if (s == zero) {
            m_data.clear();
        }
        else {

            for (DIMN i=0; i<m_dimension; ++i) {
                m_data[i] *= s;
            }
        }
        return *this;
    }

    /// Inplace rational division
    dense_vector& operator/=(const RATIONAL s)
    {
        // Instead of using the /= operator, us the *=
        // The compiler probably applies this optimisation
        // automatically, but let's be safe.
        SCALAR val = one / s;
        for (DIMN i=0; i<m_dimension; ++i) {
            m_data[i] *= val;
        }
        return *this;
    }

    /// Scalar multiplication
    __DECLARE_BINARY_OPERATOR(dense_vector, *, *=, SCALAR);

    /// Rational division
    __DECLARE_BINARY_OPERATOR(dense_vector, /, /=, RATIONAL);

    /// Inplace addition
    dense_vector& operator+=(const dense_vector& rhs)
    {
        if (empty())
            return (*this = rhs);

        if (rhs.empty())
            return *this;

        if (rhs.m_dimension > m_dimension) {
            resize_to_dimension(rhs.m_dimension);
        }

        for (DIMN i=0; i<rhs.m_dimension; ++i) {
            m_data[i] += rhs.m_data[i];
        }

        return *this;
    }

    /// Inplace subtraction
    dense_vector& operator-=(const dense_vector& rhs)
    {
        if (rhs.empty())
            return *this;

        if (rhs.m_dimension > m_dimension) {
            resize_to_dimension(rhs.m_dimension);
        }

        for (DIMN i=0; i<rhs.m_dimension; ++i) {
            m_data[i] -= rhs.m_data[i];
        }

        return *this;
    }

    /// Addition
    __DECLARE_BINARY_OPERATOR(dense_vector, +, +=, dense_vector);

    /// Subtraction
    __DECLARE_BINARY_OPERATOR(dense_vector, -, -=, dense_vector);

    /// Inplace coordinatewise min
    dense_vector& operator&=(const dense_vector& rhs)
    {
        DIMN mid = std::min(m_dimension, rhs.m_dimension);

        for (DIMN i=0; i<mid; ++i) {
            m_data[i] = std::min(m_data[i], rhs.m_data[i]);
        }

        if (m_dimension < rhs.m_dimension) {
            resize_to_dimension(rhs.m_dimension);

            for (DIMN i=mid; i<rhs.m_dimension; ++i) {
                m_data[i] = std::min(zero, rhs.m_data[i]);
            }
        }
        else if (m_dimension > rhs.m_dimension) {

            for (DIMN i=mid; i<m_dimension; ++i) {
                m_data[i] = std::min(zero, rhs.m_data[i]);
            }
        }
        return *this;
    }

    /// Inplace coordinatewise max
    dense_vector& operator|=(const dense_vector& rhs)
    {
        DIMN mid = std::min(m_dimension, rhs.m_dimension);

        for (DIMN i=0; i<mid; ++i) {
            m_data[i] = std::max(m_data[i], rhs.m_data[i]);
        }

        if (m_dimension < rhs.m_dimension) {
            resize_to_dimension(rhs.m_dimension);

            for (DIMN i=mid; i<rhs.m_dimension; ++i) {
                m_data[i] = std::max(zero, rhs.m_data[i]);
            }
        }

            for (DIMN i=mid; i<m_dimension; ++i) {
                m_data[i] = std::max(zero, m_data[i]);
            }
        }
        return *this;
    }

    /// Coordinatewise min
    __DECLARE_BINARY_OPERATOR(dense_vector, &, &=, dense_vector);

    /// Coordinatewise max
    __DECLARE_BINARY_OPERATOR(dense_vector, |, |=, dense_vector);


public:

#define DECLARE_FUSED_OP(NAME, OP1, OP2, T)                           \
    dense_vector& NAME(const KEY& rhs, const T s) {                   \
        operator[](rhs) OP1 (one OP2 s);                              \
        return *this;                                                 \
    }                                                                 \
                                                                      \
    dense_vector& NAME(const dense_vector& rhs, const T s) {          \
        if (m_dimension < rhs.m_dimension) {                          \
            resize_to_diension(rhs.m_dimension);                      \
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

    bool operator==(const dense_vector& rhs) const
    {
        DIMN mid = std::min(m_dimension, rhs.m_dimension);

        for (DIMN i=0; i<mid; ++i) {
            if (m_data[i] != rhs.m_data[i]) {
                return false;
            }
        }

        for (DIMN i=mid; i<m_dimension; ++i) {
            if (m_data[i] != zero) {
                return false;
            }
        }

        for (DIMN i=mid; i<rhs.m_dimension; ++i) {
            if (rhs.m_data[i] != zero) {
                return false;
            }
        }
        return true;
    }

    bool operator!=(const dense_vector& rhs) const
    {
        return !operator==(rhs);
    }

public:

    // Norms

    SCALAR NormL1() const
    {
        SCALAR ans(zero);
        for (DIMN i=0; i<m_dimension; ++i) {
            ans += abs(m_data[i]);
        }
        return ans;
    }

    SCALAR NormL1(const DEG deg) const
    {
        SCALAR ans(zero);


        return ans;
    }

public:

    static bool comp(typename std::pair<KEY, SCALAR> lhs,
                     typename std::pair<KEY, SCALAR> rhs)
    {
        return lhs.first < rhs.first;
    }

public:

    inline friend std::ostream& operator<<(std::ostream& os,
            const dense_vector& rhs)
    {
        std::pair<BASIS*, KEY> token;
        token.first = &dense_vector::basis;

        os << '{';
        for (DIMN i=0; i<m_dimension; ++i) {
            token.second = key_from_index(i);
            os << ' ' << m_data[i] << '(' << token << ')';
        }
        os << ' ' << '}';
        return os;
    }


};

} // namespace vectors
} // namespace alg

#endif //LIBALGEBRA_DENSE_VECTOR_H
