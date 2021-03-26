//
// Created by sam on 08/02/2021.
//


#ifndef LIBALGEBRA_DENSE_VECTOR_H
#define LIBALGEBRA_DENSE_VECTOR_H

#include <vector>
#include <utility>

#include "libalgebra/vectors/base_vector.h"
#include "libalgebra/utils/meta.h"
#include "libalgebra/basis/basis.h"
#include "libalgebra/vectors/vector.h"
#include "libalgebra/vectors/iterators.h"

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






template<typename Basis, typename Coeffs,
        typename Storage = std::vector<typename Coeffs::S> >
class dense_vector : protected base_vector<Basis, Coeffs>, dtl::requires_order<Basis>
{
    typedef Storage STORAGE;
    typedef base_vector <Basis, Coeffs> BASE_VEC;

private:
    // Data members
    STORAGE m_data;
    DIMN m_dimension;
    DEG m_degree;

public:

    // Type definitions
    typedef Basis BASIS;
    typedef Coeffs COEFFS;
    typedef typename BASIS::KEY KEY;
    typedef typename COEFFS::S SCALAR;
    typedef typename COEFFS::Q RATIONAL;

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
    dense_vector() : m_data(), m_dimension(0), m_degree(0)
    {}

    /// Copy constructor
    dense_vector(const dense_vector &v)
            : m_data(v.m_data),
              m_dimension(v.m_dimension),
              m_degree(v.m_degree)
    {}

    /// Unidimensional constructor
    explicit dense_vector(const KEY &k, const SCALAR &s = one)
            : m_dimension(0), m_degree(0)
    {
        m_data = STORAGE();
        DIMN idx = resize_for_key(k, degree_tag);
        assert (m_dimension == m_data.size());
        assert (m_data.size() > idx);
        m_data[idx] = s;
        set_degree(degree_tag);
    }

private:

    template<DEG D>
    DIMN resize_for_key(const KEY &key, alg::basis::with_degree<D>)
    {
        DEG d = basis.degree(key);
        resize_to_degree(d);
        return key_to_index(key);
    }

    DIMN resize_for_key(const KEY &key, alg::basis::without_degree)
    {
        DIMN idx = key_to_index(key);
        resize_to_dimension(idx + 1);
        return idx;
    }

    template<DEG D>
    DIMN adjust_dimension(const DIMN dim, alg::basis::with_degree<D>) const
    {
        if (dim >= max_dimension(degree_tag)) {
            return max_dimension(degree_tag);
        }

        DEG d = index_to_degree(dim);
        if (dim == start_of_degree(d)) {
            return dim;
        }
        assert (d <= D);
        return start_of_degree(d + 1);
    }

    DIMN adjust_dimension(const DIMN dim, alg::basis::without_degree) const
    {
        return std::min(max_dimension(degree_tag), dim);
    }

    template<DEG D>
    void set_degree(alg::basis::with_degree<D>)
    {
        if (m_dimension == 0) {
            m_degree = 0;
        } else {
            m_degree = index_to_degree(m_dimension - 1);
        }
        assert(m_degree <= D);
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
        assert(new_dim <= max_dimension(degree_tag));
        m_data.resize(new_dim, zero);
        m_dimension = new_dim;
        assert (m_dimension == m_data.size());
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
        assert (m_dimension == m_data.size());
        m_degree = target_deg;
    }

    DIMN next_resize_size() const
    {
        return adjust_dimension(dimension()+1, degree_tag);
    }


public:

    class iterator_item
    {
        friend class iterators::vector_iterator<iterator_item>;

        friend class dense_vector;

    public:

        iterator_item() : m_vector(NULL), m_iterator()
        {}

        iterator_item(const iterator_item &other)
                : m_vector(other.m_vector), m_iterator(other.m_iterator)
        {}

        iterator_item(dense_vector &vect, typename STORAGE::iterator it)
                : m_vector(&vect), m_iterator(it)
        {}

        typedef KEY key_type;
        typedef SCALAR &value_type;

        key_type key()
        {
            assert(index() < m_vector->m_data.size());
            return index_to_key(index());
        }

        value_type value()
        {
            assert(index() < m_vector->m_data.size());
            return *m_iterator;
        }

    private:

        dense_vector *m_vector;
        typename STORAGE::iterator m_iterator;

    private:

        bool compare_iterators(const iterator_item &other) const
        {
            return (m_iterator == other.m_iterator);
        }

        void advance()
        {
            ++m_iterator;
        }

    private:
        DIMN index() const
        {
            return DIMN(std::distance(m_vector->m_data.begin(), m_iterator));
        }
    };

    class const_iterator_item
    {
        friend class iterators::vector_iterator<const_iterator_item>;

        friend class dense_vector;

    public:

        const_iterator_item() : m_vector(NULL), m_iterator()
        {}

        const_iterator_item(const const_iterator_item &other)
                : m_vector(other.m_vector), m_iterator(other.m_iterator)
        {}

        const_iterator_item(const dense_vector &vect, typename STORAGE::const_iterator it)
                : m_vector(&vect), m_iterator(it)
        {}

        typedef KEY key_type;
        typedef const SCALAR &value_type;

        key_type key()
        {
            assert(index() < m_vector->m_data.size());
            return index_to_key(index());
        }

        value_type value()
        {
            assert(index() < m_vector->m_data.size());
            return *m_iterator;
        }

    private:

        const dense_vector *m_vector;
        typename STORAGE::const_iterator m_iterator;

    private:

        bool compare_iterators(const const_iterator_item &other) const
        {
            return (m_iterator == other.m_iterator);
        }

        void advance()
        {
            ++m_iterator;
        }

    private:

        DIMN index() const
        {
            return DIMN(std::distance(m_vector->m_data.begin(), m_iterator));
        }

    };

    typedef iterators::vector_iterator <iterator_item> iterator;
    typedef iterators::vector_iterator <const_iterator_item> const_iterator;

    // iterator methods

    iterator begin()
    {
        return iterator(*this, m_data.begin());
    }

    iterator end()
    {
        return iterator(*this, m_data.end());
    }

    const_iterator begin() const
    {
        return const_iterator(*this, m_data.begin());
    }

    const_iterator end() const
    {
        return const_iterator(*this, m_data.end());
    }

    const_iterator cbegin() const
    {
        return const_iterator(*this, m_data.begin());
    }

    const_iterator cend() const
    {
        return const_iterator(*this, m_data.end());
    }

    void insert(const KEY &key, SCALAR val)
    {
        DIMN idx = key_to_index(key);
        if (idx >= m_dimension) {
            resize_to_dimension(idx + 1);
        }
        m_data[idx] = val;
    }

    std::pair<iterator, bool> insert(std::pair<const KEY, SCALAR> &arg)
    {
        DIMN idx = key_to_index(arg.first);
        std::pair<iterator, bool> rv;

        if (idx < m_dimension && m_data[idx] != zero) {
            rv.first = iterator(*this, m_data.begin() + idx);
            rv.second = false;
            return rv;
        } else {
            resize_to_dimension(idx + 1);
        }
        assert(idx < m_data.size());
        rv.first = iterator(*this, m_data.begin() + idx);
        rv.second = false;

        rv.second = true;
        m_data[idx] = arg.second;
        return rv;
    }

    template<typename InputIterator>
    void insert(InputIterator begin, InputIterator end)
    {
        typedef std::vector<std::pair<DIMN, SCALAR> > TMPVEC;
        TMPVEC tmp;
        tmp.reserve(end - begin);

        DIMN curr, max_index = 0;

        for (InputIterator it = begin; it != end; ++it) {
            curr = key_to_index(it->first);
            max_index = (max_index < curr) ? curr : max_index;
            tmp.push_back(std::pair<DIMN, SCALAR>(curr, it->second));
        }

        resize_to_dimension(max_index + 1);

        for (typename TMPVEC::const_iterator cit = tmp.begin(); cit != tmp.end(); ++cit) {
            m_data[cit->first] = cit->second;
        }

    }

    const_iterator find(const KEY &key) const
    {
        DIMN idx;
        if ((idx = key_to_index(key)) < m_dimension) {
            return const_iterator(*this, m_data.begin() + idx);
        }
        return end();
    }

    iterator find(const KEY &key)
    {
        DIMN idx;
        if ((idx = key_to_index(key)) < m_dimension) {
            return iterator(*this, m_data.begin() + idx);
        }
        return end();
    }

    iterator find(DIMN idx)
    {
        return iterator(*this, m_data.begin() + idx);
    }

    const_iterator find(DIMN idx) const
    {
        return const_iterator(*this, m_data.begin() + idx);
    }

    void erase(const KEY &key)
    {
        DIMN idx;
        if ((idx = key_to_index(key)) < m_dimension) {
            m_data[idx] = zero;
        }
    }

    void erase(iterator &it)
    {
        assert(it != end());
        m_data[it->index()] = zero;
    }


public:

    void clear()
    {
        for (DIMN i = 0; i < m_dimension; ++i) {
            m_data[i] = zero;
        }
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
    {
        assert (m_degree <= BASIS::MAX_DEGREE);
        return m_degree;
    }

    bool degree_equals(const DEG degree) const
    {
        return m_degree == degree;
    }

    SCALAR &value(const DIMN dim)
    {
        assert (dim < m_data.size());
        return m_data[dim];
    }

    const SCALAR &value(const DIMN dim) const
    {
        assert (dim < m_data.size());
        return m_data[dim];
    }

    DIMN size() const
    {
        DIMN sz = 0;
        assert (m_dimension == m_data.size());
        for (DIMN i = 0; i < m_dimension; ++i) {
            sz += ((m_data[i] != zero) ? 1 : 0);
        }
        return sz;
    }

    bool empty() const
    {
        if (m_data.empty()) {
            return true;
        }
        assert (m_dimension == m_data.size());
        for (DIMN i = 0; i < m_dimension; ++i) {
            if (m_data[i] != zero) {
                return false;
            }
        }
        return true;
    }

protected:

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

    template<DEG D>
    static DIMN max_dimension(alg::basis::with_degree<D>)
    {
        return basis.start_of_degree(degree_tag.max_degree + 1);
    }

    static DIMN max_dimension(alg::basis::without_degree)
    {
        return basis.max_dimension();
    }

public:

    // Element access (map API)

    const SCALAR &operator[](const KEY &k) const
    {
        DIMN idx;
        assert (m_dimension == m_data.size());
        if ((idx = key_to_index(k)) < m_dimension) {
            return value(idx);
        }
        return zero;
    }

    SCALAR &operator[](const KEY &k)
    {
        DIMN idx = key_to_index(k);
        assert (m_dimension == m_data.size());
        if ((idx < m_dimension)) {
            return value(idx);
        }

        resize_to_dimension(idx + 1);
        return value(idx);
    }


public:

    // Arithmetic operation

    /// Unary negation
    dense_vector operator-(void) const
    {
        dense_vector result;
        result.resize_to_dimension(m_dimension);
        assert (m_dimension == m_data.size());
        assert (result.m_dimension == result.m_data.size());
        for (DIMN i = 0; i < m_dimension; ++i) {
            result.m_data[i] = -m_data[i];
        }

        return result;
    }

    /// Inplace scalar multiplication
    dense_vector &operator*=(const SCALAR& s)
    {

        assert (m_dimension == m_data.size());
        for (DIMN i = 0; i < m_dimension; ++i) {
            m_data[i] *= s;
        }

        return *this;
    }

    /// Inplace rational division
    dense_vector &operator/=(const RATIONAL& s)
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
    __DECLARE_BINARY_OPERATOR(dense_vector, *, *=, SCALAR);

    /// Rational division
    __DECLARE_BINARY_OPERATOR(dense_vector,
                              /, /=, RATIONAL);

    /// Inplace addition
    dense_vector &operator+=(const dense_vector &rhs)
    {
        if (empty()) {
            *this = rhs;
            return *this;
        }

        if (rhs.empty()) {
            return *this;
        }

        if (rhs.m_dimension > m_dimension) {
            resize_to_dimension(rhs.m_dimension);
        }

        assert (m_dimension == m_data.size());
        assert (m_data.size() >= rhs.m_data.size());

        for (DIMN i = 0; i < rhs.m_dimension; ++i) {
            m_data[i] += rhs.m_data[i];
        }

        return *this;
    }

    /// Inplace subtraction
    dense_vector &operator-=(const dense_vector &rhs)
    {
        if (rhs.empty()) {
            return *this;
        }

        if (rhs.m_dimension > m_dimension) {
            resize_to_dimension(rhs.m_dimension);
        }

        assert (m_data.size() >= rhs.m_data.size());

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
                m_data[i] = std::min(zero, m_data[i]);
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
    __DECLARE_BINARY_OPERATOR(dense_vector, &, &=, dense_vector);

    /// Coordinatewise max
    __DECLARE_BINARY_OPERATOR(dense_vector,
                              |, |=, dense_vector);


public:


    // Fused add/sub scalar multiplication/division
    DECLARE_FUSED_OP(add_scal_prod, +=, *, SCALAR);

    DECLARE_FUSED_OP(sub_scal_prod, -=, *, SCALAR);

    DECLARE_FUSED_OP(add_scal_div, +=, /, RATIONAL);

    DECLARE_FUSED_OP(sub_scal_div, -=, /, RATIONAL);



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

protected:

    std::pair<DIMN, bool> equal_to_min(const dense_vector &rhs) const
    {
        DIMN mid = std::min(m_dimension, rhs.m_dimension);

        for (DIMN i = 0; i < mid; ++i) {
            if (m_data[i] != rhs.m_data[i]) {
                return std::pair<DIMN, bool>(mid, false);
            }
        }
        return std::pair<DIMN, bool>(mid, true);
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
        return ans;
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

    static bool comp(std::pair<KEY, SCALAR> lhs,
                     std::pair<KEY, SCALAR> rhs)
    {
        return lhs.first < rhs.first;
    }

public:

    inline friend std::ostream &operator<<(std::ostream &os,
                                           const dense_vector &rhs)
    {
        std::pair<BASIS *, KEY> token;
        token.first = &dense_vector::basis;
        assert (rhs.m_dimension == rhs.m_data.size());
        os << '{';
        for (DIMN i = 0; i < rhs.m_dimension; ++i) {
            if (zero != rhs.m_data[i]) {
                token.second = index_to_key(i);
                os << ' ' << rhs.m_data[i] << '(' << token << ')';
            }
        }
        os << ' ' << '}';
        return os;
    }

public:

    // Transform methods

    template<typename Vector, typename KeyTransform>
    void triangular_buffered_apply_binary_transform(
            Vector &result,
            const dense_vector &rhs,
            KeyTransform key_transform,
            const DEG max_depth
    ) const
    {
        if (empty() || rhs.empty()) {
            return;
        }

        const IDEG max_degree = static_cast<IDEG>(std::min(max_depth, m_degree + rhs.m_degree));
        dense_vector &d_result = dtl::vector_base_access::convert(result);
        d_result.resize_to_degree(static_cast<DEG>(max_degree));

        IDEG lhs_deg_min, lhs_deg_max, rhs_deg;

        for (IDEG out_deg = max_degree; out_deg >= 0; --out_deg) {
            lhs_deg_min = std::max(IDEG(0), out_deg - static_cast<IDEG>(rhs.m_degree));
            lhs_deg_max = std::min(out_deg, static_cast<IDEG>(m_degree));
            for (IDEG lhs_deg = lhs_deg_max; lhs_deg >= lhs_deg_min; --lhs_deg) {
                rhs_deg = out_deg - lhs_deg;


                assert (start_of_degree(lhs_deg+1) <= m_data.size());
                assert (start_of_degree(rhs_deg+1) <= rhs.m_data.size());
                assert (d_result.m_data.size() >= start_of_degree(out_deg+1));

                for (DIMN i = start_of_degree(lhs_deg);
                     i < start_of_degree(lhs_deg + 1);
                     ++i) {
                    for (DIMN j = start_of_degree(rhs_deg);
                         j < start_of_degree(rhs_deg + 1);
                         ++j) {
                        if (m_data[i] != zero && rhs.m_data[j] != zero) {
                            key_transform(
                                    result,
                                    index_to_key(i),
                                    m_data[i],
                                    index_to_key(j),
                                    rhs.m_data[j]
                            );
                        }
                    }
                }

            }
        }

    }

    template<typename Vector, typename KeyTransform, typename IndexTransform>
    void triangular_buffered_apply_binary_transform(
            Vector &result,
            const dense_vector &rhs,
            KeyTransform key_transform,
            IndexTransform index_transform,
            const DEG max_depth
    ) const
    {
        if (empty() || rhs.empty()) {
            return;
        }

        dense_vector &d_result = dtl::vector_base_access::convert(result);

        const DEG max_degree = std::min(max_depth, m_degree + rhs.m_degree);
        d_result.resize_to_degree(max_degree);

        assert (d_result.m_data.size() == start_of_degree(max_degree+1) );


        IDEG lhs_deg_min, lhs_deg_max, rhs_deg;

        for (IDEG out_deg = max_degree; out_deg >= 0; --out_deg) {
            lhs_deg_min = std::max(IDEG(0), out_deg - static_cast<IDEG>(rhs.degree()));
            lhs_deg_max = std::min(out_deg, static_cast<IDEG>(degree()));
            for (IDEG lhs_deg = lhs_deg_max; lhs_deg >= lhs_deg_min; --lhs_deg) {
                rhs_deg = out_deg - lhs_deg;

                DIMN lh_deg_start = start_of_degree(static_cast<DEG>(lhs_deg));
                DIMN rh_deg_start = start_of_degree(static_cast<DEG>(rhs_deg));

                assert (start_of_degree(lhs_deg+1) <= m_data.size());
                assert (start_of_degree(rhs_deg+1) <= rhs.m_data.size());
                assert (d_result.m_data.size() >= start_of_degree(out_deg+1));

                index_transform(
                        &d_result.m_data[start_of_degree(out_deg)],
                        &m_data[start_of_degree(lhs_deg)],
                        &rhs.m_data[start_of_degree(rhs_deg)],
                        start_of_degree(lhs_deg + 1) - lh_deg_start,
                        start_of_degree(rhs_deg + 1) - rh_deg_start
                );
            }
        }
    }

    template <typename KeyTransform>
    void triangular_unbuffered_apply_binary_transform(
            const dense_vector& rhs,
            KeyTransform key_transform,
            const DEG max_depth
    )
    {
        dense_vector result;
        triangular_buffered_apply_binary_transform(result, rhs, key_transform, max_depth);
        swap(result);
    }

    template <typename KeyTransform, typename IndexTransform>
    void triangular_unbuffered_apply_binary_transform(
            const dense_vector& rhs,
            KeyTransform key_transform,
            IndexTransform index_transform,
            const DEG max_depth
            )
    {
        if (empty() || rhs.empty()) {
            clear();
            return;
        }
        /*
         * Ok, there are some details to work through here.
         *
         * If the basis admits a degree 0, then it must have exactly one
         * dimension, otherwise the modified degree must be visited several times,
         * which obviously won't work. If this single element is 0 or if the basis
         * does not admit a degree 0, then we don't need to touch it at all,
         * which might lead to some speed improvements.
         *
         */

        const DIMN degree_difference_1_0 = start_of_degree(1) - start_of_degree(0);

        if (degree_difference_1_0 > 1) {
            // If degree 0 has more than 1 dimension, we have to use buffering.
            // I don't expect this happens (ever). Except perhaps some contrived
            // test examples.
            dense_vector result;
            triangular_buffered_apply_binary_transform(result, rhs, key_transform, index_transform, max_depth);
            swap(result);
            return;
        }

        const DEG old_lhs_deg = degree();
        const DEG max_degree = std::min(max_depth, m_degree + rhs.m_degree);

        resize_to_degree(max_degree);
        assert(m_data.size() == start_of_degree(max_degree+1));

        if (max_degree == 0) {
            index_transform(
                    &m_data[start_of_degree(0)],
                    &m_data[start_of_degree(0)],
                    &rhs.m_data[start_of_degree(0)],
                    degree_difference_1_0,
                    degree_difference_1_0,
                    true
            );
            return;
        }

        IDEG lhs_deg_min, lhs_deg_max, rhs_deg, offset=0;
        bool assign = true;

        for (IDEG out_deg = max_degree; out_deg >= 0; --out_deg) {
            lhs_deg_min = std::max(IDEG(0), out_deg - static_cast<IDEG>(rhs.degree()));
            assign = true;

            if (degree_difference_1_0 == 0) {
                // Basis does not admit a degree 0.
                offset = 1;
            } else if (degree_difference_1_0 == 1) {
                if (rhs.m_data[0] == one && out_deg <= old_lhs_deg) {
                    assign = false;
                    offset = 1;
                }
            }

            lhs_deg_max = std::min(out_deg, static_cast<IDEG>(old_lhs_deg));

            for (IDEG lhs_deg = lhs_deg_max - offset; lhs_deg >= lhs_deg_min; --lhs_deg) {
                rhs_deg = out_deg - lhs_deg;
                DIMN lh_deg_start = start_of_degree(static_cast<DEG>(lhs_deg));
                DIMN rh_deg_start = start_of_degree(static_cast<DEG>(rhs_deg));

                assert (start_of_degree(lhs_deg+1) <= m_data.size());
                assert (start_of_degree(rhs_deg+1) <= rhs.m_data.size());
                assert (m_data.size() >= start_of_degree(out_deg+1));

                index_transform(
                        &m_data[start_of_degree(out_deg)],
                        &m_data[start_of_degree(lhs_deg)],
                        &rhs.m_data[start_of_degree(rhs_deg)],
                        start_of_degree(lhs_deg + 1) - lh_deg_start,
                        start_of_degree(rhs_deg + 1) - rh_deg_start,
                        assign
                );

                assign = false;
            }


        }
    }


    template<typename Vector, typename KeyTransform>
    void square_buffered_apply_binary_transform(
            Vector &result,
            const dense_vector &rhs,
            KeyTransform key_transform
    ) const
    {
        if (empty() || rhs.empty()) {
            return;
        }

        dense_vector &d_result = dtl::vector_base_access::convert(result);
        d_result.resize_to_dimension(std::max(m_dimension, rhs.m_dimension));

        assert (d_result.m_data.size() >= std::max(m_dimension, rhs.m_dimension));

        for (DIMN i = 0; i < m_dimension; ++i) {
            for (DIMN j = 0; j < rhs.m_dimension; ++j) {
                key_transform(
                        result,
                        index_to_key(i),
                        m_data[i],
                        index_to_key(j),
                        rhs.m_data[j]
                );
            }
        }
    }

    template<typename Vector, typename KeyTransform, typename IndexTransform>
    void square_buffered_apply_binary_transform(
            Vector &result,
            const dense_vector &rhs,
            KeyTransform /*key_transform*/,
            IndexTransform index_transform
    ) const
    {
        if (empty() || rhs.empty()) {
            return;
        }

        dense_vector &d_result = dtl::vector_base_access::convert(result);
        d_result.resize_to_dimension(std::max(m_dimension, rhs.m_dimension));

        assert (d_result.m_data.size() >= std::max(m_dimension, rhs.m_dimension));

        index_transform(
                &d_result.m_data[0],
                &m_data[0],
                &rhs.m_data[0],
                m_dimension,
                rhs.m_dimension
        );
    }

public:

    template <typename Transform>
    void buffered_apply_unary_transform(
            dense_vector& result,
            Transform transform,
            const DEG max_deg
    ) const
    {
        if (empty()) {
            return;
        }

        result.resize_to_dimension(transform.dense_resize(dimension()));
        typename Transform::index_transform it(transform.get_index_transform());

        it(
                &result.m_data[0],
                result.dimension(),
                &m_data[0],
                dimension(),
                max_deg);
    }


    template <typename Transform>
    void buffered_apply_unary_transform(
            dense_vector& result,
            Transform transform
    ) const
    {
        if (empty()) {
            return;
        }

        result.resize_to_dimension(transform.dense_resize(dimension()));
        typename Transform::index_transform it(transform.get_index_transform());

        it(
                &result.m_data[0],
                result.dimension(),
                &m_data[0],
                dimension());
    }





};


#undef DECLARE_FUSED_OP


} // namespace vectors
} // namespace alg

#endif //LIBALGEBRA_DENSE_VECTOR_H
