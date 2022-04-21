//
// Created by sam on 08/02/2021.
//

#ifndef LIBALGEBRA_DENSE_VECTOR_H
#define LIBALGEBRA_DENSE_VECTOR_H

#include <utility>
#include <vector>

#include "libalgebra/basis/basis.h"
#include "libalgebra/utils/meta.h"
#include "libalgebra/vectors/base_vector.h"
#include "libalgebra/vectors/dense_storage.h"
#include "libalgebra/vectors/iterators.h"
#include "libalgebra/vectors/vector.h"

namespace alg {
namespace vectors {

namespace dtl {

template<typename Basis>
struct requires_order {
    typedef typename alg::basis::basis_traits<Basis>::ordering_tag::order key_ordering;
};

}// namespace dtl

/**
 * @brief Vector type where data is stored in contiguous memory.
 *
 * A dense vector stores coefficients in a contiguous block of memory, where the association with
 * basis element is formed by taking the key corresponding to the index of the coefficient in the
 * basis total ordering. Consequently, a dense vector is valid if and only if Basis is totally
 * ordered.
 *
 * @tparam Basis
 * @tparam Coeffs
 */
template<typename Basis, typename Coeffs>
class dense_vector : protected base_vector<Basis, Coeffs>, dtl::requires_order<Basis>
{
    typedef dense_storage<typename Coeffs::S> STORAGE;
    typedef base_vector<Basis, Coeffs> BASE_VEC;

    friend class dtl::data_access_base<dense_vector>;

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
    using BASE_VEC::mone;
    using BASE_VEC::one;
    using BASE_VEC::zero;

    using BASE_VEC::degree_tag;

public:
    // Constructors
    /// Default constructor
    dense_vector()
        : m_data(), m_dimension(0), m_degree(0)
    {}

    /// Copy constructor
    dense_vector(const dense_vector& v)
        : m_data(v.m_data), m_dimension(v.m_dimension), m_degree(v.m_degree)
    {}

    /// Move constructor
    dense_vector(dense_vector&& other) noexcept
        : m_data(std::move(other.m_data)), m_dimension(other.m_dimension), m_degree(other.m_degree)
    {}

    /// Unidimensional constructor
    explicit dense_vector(const KEY& k, const SCALAR& s = one)
        : m_data(), m_dimension(0), m_degree(0)
    {
        DIMN idx = resize_for_key(k, degree_tag);
        assert(m_dimension == m_data.size());
        assert(m_data.size() > idx);
        m_data[idx] = s;
        set_degree(degree_tag);
    }

    /**
     * @brief Construct from pointer to data
     * @param begin
     * @param end
     */
    dense_vector(SCALAR const* begin, SCALAR const* end)
        : m_data(begin, end),
          m_dimension(m_data.size()),
          m_degree(0)
    {
        if (m_data.size() != adjust_dimension(m_data.size(), degree_tag)) {
            resize_to_dimension(m_data.size());
        }
        else {
            set_degree(degree_tag);
        }
    }

    /**
     * @brief Create a new dense vector from data with an initial offset
     *
     * Construct a new vector in which the first offset elements are zero,
     * and fill the next (end - begin) values by copying from the range
     * [begin, end).
     *
     * @param offset Number of terms in the basis before the data starts
     * @param begin start of data range to construct
     * @param end first past last of data range to
     */
    dense_vector(DIMN offset, SCALAR const* begin, SCALAR const* end)
        : m_data(offset, begin, end),
          m_dimension(0),
          m_degree(0)
    {
        if (m_data.size() != adjust_dimension(m_data.size(), degree_tag)) {
            resize_to_dimension(m_data.size());
        }
        else {
            set_degree(degree_tag);
        }
    }

    /**
     * @brief Create a new dense vector from data with an initial offset
     *
     * Construct a new vector in which the first offset elements are zero,
     * and fill the next (end - begin) values by copying from the range
     * [begin, end).
     *
     * @param offset Number of terms in the basis before the data starts
     * @param begin start of data range to construct
     * @param end first past last of data range to
     */
    dense_vector(DIMN offset, SCALAR* begin, SCALAR* end)
        : m_data(offset, begin, end),
          m_dimension(0),
          m_degree(0)
    {
        if (m_data.size() != adjust_dimension(m_data.size(), degree_tag)) {
            resize_to_dimension(m_data.size());
        }
    }

    dense_vector& operator=(const dense_vector& other) = default;
    dense_vector& operator=(dense_vector&& other) noexcept = default;

private:
    template<DEG D>
    DIMN resize_for_key(const KEY& key, alg::basis::with_degree<D>)
    {
        DEG d = basis.degree(key);
        if (m_degree == 0 || m_degree < d) {
            resize_to_degree(d);
        }
        return key_to_index(key);
    }

    DIMN resize_for_key(const KEY& key, alg::basis::without_degree)
    {
        DIMN idx = key_to_index(key);
        if (m_data.size() == 0 || m_data.size() <= idx) {
            resize_to_dimension(idx + 1);
        }
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
        assert(d <= D);
        return start_of_degree(d + 1);
    }

    DIMN adjust_dimension(const DIMN dim, alg::basis::without_degree) const
    {
        return std::min(max_dimension(degree_tag), dim);
    }

    template<DEG D>
    void set_degree(alg::basis::with_degree<D>)
    {
        if (dimension() == 0) {
            m_degree = 0;
        }
        else {
            m_degree = index_to_degree(dimension() - 1);
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
    void reserve_to_degree(const DEG deg)
    {
        DEG target_deg = std::min(degree_tag.max_degree, deg);
        DIMN target_dim = start_of_degree(target_deg + 1);
        m_data.reserve(target_dim);
        m_dimension = target_dim;
        m_degree = target_deg;
    }

    /// Resize to dimension
    void resize_to_dimension(const DIMN dim)
    {
        DIMN new_dim = adjust_dimension(dim, degree_tag);
        assert(new_dim <= max_dimension(degree_tag));
        m_data.resize(new_dim, zero);
        m_dimension = new_dim;
        set_degree(degree_tag);
    }

    /// Reserve to degree
    void resize_to_degree(const DEG deg)
    {
        DEG target_deg = std::min(degree_tag.max_degree, deg);
        DIMN dim = start_of_degree(target_deg + 1);
        m_data.resize(dim, zero);
        m_dimension = dim;
        m_degree = target_deg;
    }

    /// Get the next valid size of a vector
    DIMN next_resize_size() const
    {
        return adjust_dimension(dimension() + 1, degree_tag);
    }

public:
    class iterator_item
    {
        friend class iterators::vector_iterator<iterator_item>;

        friend class dense_vector;

    public:
        iterator_item()
            : m_vector(nullptr), m_iterator()
        {}

        iterator_item(dense_vector& vect, typename STORAGE::iterator it)
            : m_vector(&vect), m_iterator(it)
        {}

        typedef KEY key_type;
        typedef SCALAR& value_type;

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
        dense_vector* m_vector;
        typename STORAGE::iterator m_iterator;

    private:
        bool compare_iterators(const iterator_item& other) const
        {
            return (m_iterator == other.m_iterator);
        }

        void advance()
        {
            ++m_iterator;
        }

    public:
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
        const_iterator_item()
            : m_vector(nullptr), m_iterator()
        {}

        const_iterator_item(const dense_vector& vect, typename STORAGE::const_iterator it)
            : m_vector(&vect),
              m_iterator(it)
        {}

        typedef KEY key_type;
        typedef const SCALAR& value_type;

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
        const dense_vector* m_vector;
        typename STORAGE::const_iterator m_iterator;

    private:
        bool compare_iterators(const const_iterator_item& other) const
        {
            return (m_iterator == other.m_iterator);
        }

        void advance()
        {
            ++m_iterator;
        }

    public:
        DIMN index() const
        {
            return DIMN(std::distance(m_vector->m_data.begin(), m_iterator));
        }
    };

    typedef iterators::vector_iterator<iterator_item> iterator;
    typedef iterators::vector_iterator<const_iterator_item> const_iterator;

    // iterator methods

    /// Iterator at start of vector
    iterator begin()
    {
        return iterator(*this, m_data.begin());
    }

    /// Iterator at (one past) end of vector
    iterator end()
    {
        return iterator(*this, m_data.end());
    }

    /// Const iterator at start of vector
    const_iterator begin() const
    {
        return const_iterator(*this, m_data.begin());
    }

    /// Const iterator at (one past) end of vector
    const_iterator end() const
    {
        return const_iterator(*this, m_data.end());
    }

    /// Const iterator at start of vector
    const_iterator cbegin() const
    {
        return const_iterator(*this, m_data.begin());
    }

    /// Const iterator at (one past) end of vector
    const_iterator cend() const
    {
        return const_iterator(*this, m_data.end());
    }

    /**
     * @brief Insert a new value at key
     * @param key key position at which to insert new value
     * @param val new value to insert
     */
    void insert(const KEY& key, SCALAR val)
    {
        DIMN idx = key_to_index(key);
        if (idx >= dimension()) {
            resize_to_dimension(idx + 1);
        }
        m_data[idx] = val;
    }

    /**
     * @brief Insert a new value by providing a key-value pair
     * @param arg key value pair with position and value to insert into the vector
     * @return pair of iterator to place where value was inserted and bool indicating whether a value was inserted.
     *        Currently this is always true.
     */
    std::pair<iterator, bool> insert(std::pair<const KEY, SCALAR>& arg)
    {
        DIMN idx = key_to_index(arg.first);
        std::pair<iterator, bool> rv;

        if (idx < dimension() && m_data[idx] != zero) {
            rv.first = iterator(*this, m_data.begin() + idx);
            rv.second = false;
            return rv;
        }
        else {
            resize_to_dimension(idx + 1);
        }
        assert(idx < m_data.size());
        rv.first = iterator(*this, m_data.begin() + idx);
        rv.second = false;

        rv.second = true;
        m_data[idx] = arg.second;
        return rv;
    }

    /**
     * @brief Insert new values from iterator
     * @tparam InputIterator Iterator type (any class supporting InputIterator trait)
     * with value_type is a pair of KEY and SCALAR.
     * @param begin start of range to insert
     * @param end one past end of range to insert
     */
    template<typename InputIterator>
    void insert(InputIterator begin, InputIterator end)
    {
        typedef std::vector<std::pair<DIMN, SCALAR>> TMPVEC;
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

    /**
     * @brief Get a const iterator to value corresponding to key
     * @param key Key to find value for
     * @return const iterator pointing to the value
     */
    const_iterator find(const KEY& key) const
    {
        DIMN idx;
        if ((idx = key_to_index(key)) < dimension()) {
            return const_iterator(*this, m_data.begin() + idx);
        }
        return end();
    }

    /**
     * @brief Get a iterator to value corresponding to key
     * @param key Key to find value for
     * @return iterator pointing to the value
     */
    iterator find(const KEY& key)
    {
        DIMN idx;
        if ((idx = key_to_index(key)) < dimension()) {
            return iterator(*this, m_data.begin() + idx);
        }
        return end();
    }

    /**
     * @brief Get an iterator to a specific index in the vector.
     * @param idx Index of value to find
     * @return iterator pointing to the value
     */
    iterator find_index(DIMN idx)
    {
        return iterator(*this, m_data.begin() + idx);
    }

    /**
     * @brief Get a const iterator to a specific index in the vector.
     * @param idx Index of value to find
     * @return const iterator pointing to the value
     */
    const_iterator find_index(DIMN idx) const
    {
        return const_iterator(*this, m_data.begin() + idx);
    }

    /**
     * @brief Erase (set to zero) the value corresponding to key
     * @param key key to erase
     */
    void erase(const KEY& key)
    {
        DIMN idx;
        if ((idx = key_to_index(key)) < dimension()) {
            m_data[idx] = zero;
        }
    }

    /**
     * @brief Erase the value at the iterator position
     * @param it iterator to value to erase
     */
    void erase(iterator& it)
    {
        assert(it != end());
        m_data[it->index()] = zero;
    }

public:
    /**
     * @brief Clear (set to zero) all values in the vector.
     */
    void clear()
    {
        SCALAR z = zero;
        for (DIMN i = 0; i < dimension(); ++i) {
            m_data[i] = z;
        }
    }

public:
    /// Swap operation
    void swap(dense_vector& other)
    {
        m_data.swap(other.m_data);
        std::swap(m_dimension, other.m_dimension);
        std::swap(m_degree, other.m_degree);
    }

public:
    // Access to information about the vector and coefficients

    /// get the dimension of the vector
    DIMN dimension() const
    {
        return m_data.size();
    }

    /// get the maximum degree element represented by the vector
    DEG degree() const
    {
        return m_degree;
    }

    /// Check if the degree of the vector is equal to given vector
    bool degree_equals(const DEG degree) const
    {
        return m_degree == degree;
    }

    /// Get a reference to the value at index
    SCALAR& value(const DIMN dim)
    {
        assert(dim < m_data.size());
        return m_data[dim];
    }

    /// Get a const reference to the value at index
    const SCALAR& value(const DIMN dim) const
    {
        assert(dim < m_data.size());
        return m_data[dim];
    }

    /**
     * @brief Get the number of non-zero elements in the vector
     *
     * This is a fairly slow operation since it must check every element in the vector.
     *
     * @return Number of non-zero elements.
     */
    DIMN size() const
    {
        DIMN sz = 0;
        for (DIMN i = 0; i < dimension(); ++i) {
            sz += ((m_data[i] != zero) ? 1 : 0);
        }
        return sz;
    }

    /**
     * @brief Check if the vector contains only zero elements.
     *
     * This is a fairly slow operation, especially if there are a large number of zeros.
     *
     * @return true if vector is empty (only has zeros) and false otherwise
     */
    bool empty() const
    {
        if (m_data.empty()) {
            return true;
        }

        for (DIMN i = 0; i < dimension(); ++i) {
            if (m_data[i] != zero) {
                return false;
            }
        }
        return true;
    }

protected:
    /// Index of key and key of index
    static KEY index_to_key(const DIMN index)
    {
        return basis.index_to_key(index);
    }

    /// Get the index corresponding to a key
    static DIMN key_to_index(const KEY& key)
    {
        return basis.key_to_index(key);
    }

    /// Get the degree of the key at an index
    static DEG index_to_degree(const DIMN idx)
    {
        return basis.degree(index_to_key(idx));
    }

    /// Get the index at which the elements of degree start
    static DIMN start_of_degree(const DEG deg)
    {
        return basis.start_of_degree(deg);
    }

    /// Get he maximum feasible dimension for bases with degree
    template<DEG D>
    static DIMN max_dimension(alg::basis::with_degree<D>)
    {
        return basis.start_of_degree(degree_tag.max_degree + 1);
    }

    /// Get the maximum feasible dimension for a basis without degree
    static DIMN max_dimension(alg::basis::without_degree)
    {
        return basis.max_dimension();
    }

public:
    // Element access (map API)

    const SCALAR& operator[](const KEY& k) const
    {
        DIMN idx;
        assert(m_dimension == m_data.size());
        if ((idx = key_to_index(k)) < dimension()) {
            return value(idx);
        }
        return zero;
    }

    SCALAR& operator[](const KEY& k)
    {
        DIMN idx = key_to_index(k);

        if ((idx < dimension())) {
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
        result.resize_to_dimension(dimension());
        assert(m_dimension == m_data.size());
        assert(result.m_dimension == result.m_data.size());
        for (DIMN i = 0; i < dimension(); ++i) {
            result.m_data[i] = -m_data[i];
        }

        return result;
    }

    /// Inplace scalar multiplication
    dense_vector& operator*=(const SCALAR& s)
    {

        assert(m_dimension == m_data.size());
        for (DIMN i = 0; i < dimension(); ++i) {
            m_data[i] *= s;
        }

        return *this;
    }

    /// Inplace rational division
    dense_vector& operator/=(const RATIONAL& s)
    {
        // Instead of using the /= operator, us the *=
        // The compiler probably applies this optimisation
        // automatically, but let's be safe.
        SCALAR val = one / s;
        for (DIMN i = 0; i < dimension(); ++i) {
            m_data[i] *= val;
        }
        return *this;
    }

    /// Scalar multiplication
    dense_vector operator*(const SCALAR& rhs) const
    {
        dense_vector result(*this);
        return result *= rhs;
    };

    /// Rational division
    dense_vector operator/(const RATIONAL& rhs) const
    {
        dense_vector result(*this);
        return result /= rhs;
    };

    /// Inplace addition
    dense_vector& operator+=(const dense_vector& rhs)
    {
        if (dimension() == 0) {
            *this = rhs;
            return *this;
        }

        if (rhs.dimension() == 0) {
            return *this;
        }

        if (rhs.dimension() > dimension()) {
            resize_to_dimension(rhs.dimension());
            //m_data.copy_extend(rhs.m_data.begin() + mid_dim, rhs.m_data.end());
            assert(dimension() == rhs.dimension());
        }

        SCALAR* lh_ptr = m_data.begin();
        SCALAR const* rh_ptr = rhs.m_data.cbegin();

        //for (DIMN i = 0; i < mid_dim; ++i) {
        //   lh_ptr[i] += rh_ptr[i];
        //}
        for (DIMN i = 0; i < rhs.dimension(); ++i) {
            lh_ptr[i] += rh_ptr[i];
        }

        return *this;
    }

    /// Inplace subtraction
    dense_vector& operator-=(const dense_vector& rhs)
    {
        if (rhs.empty()) {
            return *this;
        }

        if (rhs.m_dimension > dimension()) {
            resize_to_dimension(rhs.dimension());
        }

        for (DIMN i = 0; i < rhs.dimension(); ++i) {
            m_data[i] -= rhs.m_data[i];
        }

        return *this;
    }

    /// Addition
    dense_vector operator+(const dense_vector& rhs) const
    {
        dense_vector result(*this);
        return result += rhs;
    };

    /// Subtraction
    dense_vector operator-(const dense_vector& rhs) const
    {
        dense_vector result(*this);
        return result -= rhs;
    };

    /// Inplace coordinatewise min
    dense_vector& operator&=(const dense_vector& rhs)
    {
        DIMN mid = std::min(dimension(), rhs.dimension());

        for (DIMN i = 0; i < mid; ++i) {
            m_data[i] = std::min(m_data[i], rhs.m_data[i]);
        }

        if (dimension() < rhs.dimension()) {
            resize_to_dimension(rhs.m_dimension);

            for (DIMN i = mid; i < rhs.dimension(); ++i) {
                m_data[i] = std::min(zero, rhs.m_data[i]);
            }
        }
        else if (dimension() > rhs.dimension()) {

            for (DIMN i = mid; i < dimension(); ++i) {
                m_data[i] = std::min(zero, m_data[i]);
            }
        }
        return *this;
    }

    /// Inplace coordinatewise max
    dense_vector& operator|=(const dense_vector& rhs)
    {
        DIMN mid = std::min(dimension(), rhs.dimension());

        for (DIMN i = 0; i < mid; ++i) {
            m_data[i] = std::max(m_data[i], rhs.m_data[i]);
        }

        if (dimension() < rhs.dimension()) {
            resize_to_dimension(rhs.dimension());

            for (DIMN i = mid; i < rhs.dimension(); ++i) {
                m_data[i] = std::max(zero, rhs.m_data[i]);
            }
        }

        for (DIMN i = mid; i < dimension(); ++i) {
            m_data[i] = std::max(zero, m_data[i]);
        }

        return *this;
    }

    /// Coordinatewise min
    dense_vector operator&(const dense_vector& rhs) const
    {
        dense_vector result(*this);
        return result &= rhs;
    };

    /// Coordinatewise max
    dense_vector operator|(const dense_vector& rhs) const
    {
        dense_vector result(*this);
        return result |= rhs;
    };

public:
    // Fused add/sub scalar multiplication/division
    dense_vector& add_scal_prod(const KEY& rhs, const SCALAR s)
    {
        operator[](rhs) += (one * s);
        return *this;
    }

    dense_vector& add_scal_prod(const dense_vector& rhs, const SCALAR s)
    {
        if (rhs.dimension() == 0) {
            return *this;
        }

        if (dimension() < rhs.dimension()) {
            resize_to_dimension(rhs.dimension());
        }
        SCALAR* lh_ptr = &m_data[0];
        SCALAR const* rh_ptr = &rhs.m_data[0];
        for (DIMN i = 0; i < rhs.dimension(); ++i) {
            lh_ptr[i] += (rh_ptr[i] * s);
        }
        return *this;
    };

    dense_vector& sub_scal_prod(const KEY& rhs, const SCALAR s)
    {
        operator[](rhs) -= (one * s);
        return *this;
    }

    dense_vector& sub_scal_prod(const dense_vector& rhs, const SCALAR s)
    {
        if (rhs.dimension() == 0) {
            return *this;
        }
        if (dimension() < rhs.dimension()) {
            resize_to_dimension(rhs.dimension());
        }
        SCALAR* lh_ptr = &m_data[0];
        SCALAR const* rh_ptr = &rhs.m_data[0];
        for (DIMN i = 0; i < rhs.dimension(); ++i) {
            lh_ptr[i] -= (rh_ptr[i] * s);
        }
        return *this;
    };

    dense_vector& add_scal_div(const KEY& rhs, const RATIONAL s)
    {
        operator[](rhs) += (one / s);
        return *this;
    }

    dense_vector& add_scal_div(const dense_vector& rhs, const RATIONAL s)
    {
        if (rhs.dimension() == 0) {
            return *this;
        }
        if (dimension() < rhs.dimension()) {
            resize_to_dimension(rhs.dimension());
        }
        SCALAR* lh_ptr = &m_data[0];
        SCALAR const* rh_ptr = &rhs.m_data[0];
        for (DIMN i = 0; i < rhs.dimension(); ++i) {
            lh_ptr[i] += (rh_ptr[i] / s);
        }
        return *this;
    };

    dense_vector& sub_scal_div(const KEY& rhs, const RATIONAL s)
    {
        operator[](rhs) -= (one / s);
        return *this;
    }

    dense_vector& sub_scal_div(const dense_vector& rhs, const RATIONAL s)
    {
        if (rhs.dimension() == 0) {
            return *this;
        }
        if (dimension() < rhs.dimension()) {
            resize_to_dimension(rhs.dimension());
        }
        SCALAR* lh_ptr = &m_data[0];
        SCALAR const* rh_ptr = &rhs.m_data[0];
        for (DIMN i = 0; i < rhs.dimension(); ++i) {
            lh_ptr[i] -= (rh_ptr[i] / s);
        }
        return *this;
    };

public:
    bool operator==(const dense_vector& rhs) const
    {
        DIMN mid = std::min(dimension(), rhs.dimension());

        for (DIMN i = 0; i < mid; ++i) {
            if (m_data[i] != rhs.m_data[i]) {
                return false;
            }
        }

        for (DIMN i = mid; i < dimension(); ++i) {
            if (m_data[i] != zero) {
                return false;
            }
        }

        for (DIMN i = mid; i < rhs.dimension(); ++i) {
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

protected:
    std::pair<DIMN, bool> equal_to_min(const dense_vector& rhs) const
    {
        DIMN mid = std::min(dimension(), rhs.dimension());

        for (DIMN i = 0; i < mid; ++i) {
            if (m_data[i] != rhs.m_data[i]) {
                return std::pair<DIMN, bool>(mid, false);
            }
        }
        return std::pair<DIMN, bool>(mid, true);
    }

public:
    // Norms

    /**
     * @brief Compute the L1 norm of the vector
     *
     * The L1 norm here is the sum of absolute values of members of the vector.
     *
     * @return the L1 norm
     */
    SCALAR NormL1() const
    {
        SCALAR ans(zero);
        for (DIMN i = 0; i < dimension(); ++i) {
            ans += abs(m_data[i]);
        }
        return ans;
    }

    /**
     * @brief Compute the L1 norm of elements with degree at most deg
     * @param deg maximum degree of elements
     * @return the L1 norm of elements
     */
    SCALAR NormL1(const DEG deg) const
    {
        SCALAR ans(zero);
        for (DIMN i = start_of_degree(deg); i < start_of_degree(deg + 1); ++i) {
            ans += abs(m_data[i]);
        }
        return ans;
    }

    /**
     * @brief Compute the L-infinity norm of the the vector
     *
     * The L-infinity norm is the maximum absolute value of elements in the vector.
     *
     * @return the L-infinity norm of the vector
     */
    SCALAR NormLInf() const
    {
        SCALAR ans(zero);
        for (DIMN i = 0; i < dimension(); ++i) {
            SCALAR abs_val = abs(m_data[i]);
            ans = (abs_val > ans) ? abs_val : ans;
        }
        return ans;
    }

    /**
     * @brief Compute the L-infinity norm of elements in the vector of degree at most deg
     * @param deg maximum degree of elements
     * @return the L-infinity norm of elements
     */
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
    static bool comp(std::pair<KEY, SCALAR> lhs, std::pair<KEY, SCALAR> rhs)
    {
        return lhs.first < rhs.first;
    }

protected:
    void print_members(std::ostream& os) const
    {
        std::pair<BASIS*, KEY> token;
        token.first = &dense_vector::basis;
        for (DIMN i = 0; i < dimension(); ++i) {
            if (zero != m_data[i]) {
                token.second = index_to_key(i);
                os << ' ' << m_data[i] << '(' << token << ')';
            }
        }
    }

public:
    inline friend std::ostream& operator<<(std::ostream& os, const dense_vector& rhs)
    {
        os << '{';
        rhs.print_members(os);
        os << ' ' << '}';
        return os;
    }

public:
    // Transform methods

    /**
     * @brief Apply a buffered binary transform using only key transform up to max depth
     *
     * This is applied to the vector using the degree optimisation.
     *
     * @tparam Vector Result vector type
     * @tparam KeyTransform Key transform type
     * @param result Vector in which to place the result
     * @param rhs right hand side buffer
     * @param key_transform transform to apply
     * @param max_depth maximum depth of elements to compute
     */
    template<typename Vector, typename KeyTransform>
    void
    triangular_buffered_apply_binary_transform(Vector& result, const dense_vector& rhs, KeyTransform key_transform,
                                               const DEG max_depth) const
    {
        if (empty() || rhs.empty()) {
            return;
        }

        const IDEG max_degree = static_cast<IDEG>(std::min(max_depth, m_degree + rhs.m_degree));
        dense_vector& d_result = dtl::vector_base_access::convert(result);
        d_result.resize_to_degree(static_cast<DEG>(max_degree));

        IDEG lhs_deg_min, lhs_deg_max, rhs_deg;

        for (IDEG out_deg = max_degree; out_deg >= 0; --out_deg) {
            lhs_deg_min = std::max(IDEG(0), out_deg - static_cast<IDEG>(rhs.m_degree));
            lhs_deg_max = std::min(out_deg, static_cast<IDEG>(m_degree));
            for (IDEG lhs_deg = lhs_deg_max; lhs_deg >= lhs_deg_min; --lhs_deg) {
                rhs_deg = out_deg - lhs_deg;

                assert(start_of_degree(lhs_deg + 1) <= m_data.size());
                assert(start_of_degree(rhs_deg + 1) <= rhs.m_data.size());
                assert(d_result.m_data.size() >= start_of_degree(out_deg + 1));

                for (DIMN i = start_of_degree(lhs_deg); i < start_of_degree(lhs_deg + 1); ++i) {
                    for (DIMN j = start_of_degree(rhs_deg); j < start_of_degree(rhs_deg + 1); ++j) {
                        if (m_data[i] != zero && rhs.m_data[j] != zero) {
                            key_transform(result, index_to_key(i), m_data[i], index_to_key(j), rhs.m_data[j]);
                        }
                    }
                }
            }
        }
    }

    /**
     * @brief Apply a buffered binary transform with separate transforms up to max depth
     *
     * This is applied to the vector using the degree optimisation.
     *
     * @tparam Vector Result vector type
     * @tparam KeyTransform Key transform type
     * @tparam IndexTransform Index transform type
     * @param result Vector in which to place the result
     * @param rhs Right hand side buffer
     * @param key_transform transform to apply by keys (sparse elements)
     * @param index_transform transform to apply by index (dense elements)
     * @param max_depth Maximum depth to compute
     */
    template<typename Vector, typename KeyTransform, typename IndexTransform>
    void
    triangular_buffered_apply_binary_transform(Vector& result, const dense_vector& rhs, KeyTransform,
                                               IndexTransform index_transform, const DEG max_depth) const
    {
        if (empty() || rhs.empty()) {
            return;
        }

        dense_vector& d_result = dtl::vector_base_access::convert(result);

        const DEG max_degree = std::min(max_depth, m_degree + rhs.m_degree);
        d_result.resize_to_degree(max_degree);

        assert(d_result.m_data.size() == start_of_degree(max_degree + 1));

        IDEG lhs_deg_min, lhs_deg_max, rhs_deg;

        for (IDEG out_deg = max_degree; out_deg >= 0; --out_deg) {
            lhs_deg_min = std::max(IDEG(0), out_deg - static_cast<IDEG>(rhs.degree()));
            lhs_deg_max = std::min(out_deg, static_cast<IDEG>(degree()));
            for (IDEG lhs_deg = lhs_deg_max; lhs_deg >= lhs_deg_min; --lhs_deg) {
                rhs_deg = out_deg - lhs_deg;

                DIMN lh_deg_start = start_of_degree(static_cast<DEG>(lhs_deg));
                DIMN rh_deg_start = start_of_degree(static_cast<DEG>(rhs_deg));

                assert(start_of_degree(static_cast<DEG>(lhs_deg + 1)) <= m_data.size());
                assert(start_of_degree(static_cast<DEG>(rhs_deg + 1)) <= rhs.m_data.size());
                assert(d_result.m_data.size() >= start_of_degree(static_cast<DEG>(out_deg + 1)));

                index_transform(
                        &d_result.m_data[start_of_degree(static_cast<DEG>(out_deg))],
                        &m_data[lh_deg_start],
                        &rhs.m_data[rh_deg_start],
                        start_of_degree(static_cast<DEG>(lhs_deg + 1)) - lh_deg_start,
                        start_of_degree(static_cast<DEG>(rhs_deg + 1)) - rh_deg_start);
            }
        }
    }

    /**
     * @brief Apply an unbuffered binary transform using only key transform
     *
     * @tparam KeyTransform Key transform type
     * @param rhs Right hand side buffer
     * @param key_transform Key transform
     * @param max_depth maximum depth of elements to compute
     */
    template<typename KeyTransform>
    void
    triangular_unbuffered_apply_binary_transform(const dense_vector& rhs, KeyTransform key_transform,
                                                 const DEG max_depth)
    {
        dense_vector result;
        triangular_buffered_apply_binary_transform(result, rhs, key_transform, max_depth);
        swap(result);
    }

    /**
     * @brief Apply an unbuffered binary transform with separate transforms
     * @tparam KeyTransform Key transform type
     * @tparam IndexTransform Index transform type
     * @param rhs Right hand side buffer
     * @param key_transform transform to apply by key (sparse)
     * @param index_transform transform to apply by index (dense)
     * @param max_depth Maximum depth of elements to compute
     */
    template<typename KeyTransform, typename IndexTransform>
    void
    triangular_unbuffered_apply_binary_transform(const dense_vector& rhs, KeyTransform key_transform,
                                                 IndexTransform index_transform, const DEG max_depth)
    {
        if (dimension() == 0 || rhs.dimension() == 0) {
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

        const IDEG old_lhs_deg = static_cast<IDEG>(degree());
        const DEG max_degree = std::min(max_depth, m_degree + rhs.m_degree);

        if (max_degree > m_degree) {
            //resize_to_degree(max_degree);
            reserve_to_degree(max_degree);
        }
        assert(m_data.size() >= start_of_degree(max_degree + 1));

        if (max_degree == DEG(0)) {
            index_transform(
                    &m_data[start_of_degree(0)],
                    &m_data[start_of_degree(0)],
                    &rhs.m_data[start_of_degree(0)],
                    degree_difference_1_0,
                    degree_difference_1_0,
                    true);
            return;
        }

        IDEG lhs_deg_min, lhs_deg_max, rhs_deg, offset = 0;
        bool assign, default_assign = true;

        if (degree_difference_1_0 == 1 && rhs.m_data[0] == zero) {
            offset = 1;
        }

        const IDEG max_rhs_deg = static_cast<IDEG>(rhs.degree());

        for (IDEG out_deg = static_cast<IDEG>(max_degree); out_deg >= 1; --out_deg) {
            lhs_deg_min = std::max(IDEG(0), out_deg - max_rhs_deg);
            assign = (out_deg > old_lhs_deg) || default_assign;
            lhs_deg_max = std::min(out_deg - offset, old_lhs_deg);

            for (IDEG lhs_deg = lhs_deg_max; lhs_deg >= lhs_deg_min; --lhs_deg) {
                rhs_deg = out_deg - lhs_deg;
                DIMN lh_deg_start = start_of_degree(static_cast<DEG>(lhs_deg));
                DIMN rh_deg_start = start_of_degree(static_cast<DEG>(rhs_deg));

                assert(start_of_degree(static_cast<DEG>(lhs_deg + 1)) <= m_data.size());
                assert(start_of_degree(static_cast<DEG>(rhs_deg + 1)) <= rhs.m_data.size());
                assert(m_data.size() >= start_of_degree(static_cast<DEG>(out_deg + 1)));

                index_transform(
                        &m_data[start_of_degree(static_cast<DEG>(out_deg))],
                        &m_data[lh_deg_start],
                        &rhs.m_data[rh_deg_start],
                        start_of_degree(static_cast<DEG>(lhs_deg + 1)) - lh_deg_start,
                        start_of_degree(static_cast<DEG>(rhs_deg + 1)) - rh_deg_start,
                        assign);

                assign = false;
            }
        }

        if (degree_difference_1_0 == 1) {
            index_transform(
                    &m_data[0], &m_data[0], &rhs.m_data[0], DIMN(1), DIMN(1), true);
        }
    }

    /**
     * @brief Apply buffered binary transform with no degree optimisation
     * @tparam Vector Output vector type
     * @tparam KeyTransform Key transform type
     * @param result buffer in which to place result
     * @param rhs Right hand side buffer
     * @param key_transform Transform to apply by key (sparse)
     */
    template<typename Vector, typename KeyTransform>
    void
    square_buffered_apply_binary_transform(Vector& result, const dense_vector& rhs, KeyTransform key_transform) const
    {
        if (empty() || rhs.empty()) {
            return;
        }

        dense_vector& d_result = dtl::vector_base_access::convert(result);
        d_result.resize_to_dimension(std::max(dimension(), rhs.dimension()));

        assert(d_result.m_data.size() >= std::max(dimension(), rhs.dimension()));

        for (DIMN i = 0; i < dimension(); ++i) {
            for (DIMN j = 0; j < rhs.dimension(); ++j) {
                key_transform(result, index_to_key(i), m_data[i], index_to_key(j), rhs.m_data[j]);
            }
        }
    }

    /**
     * @brief Apply buffered binary transform with separate transforms and no degree optimisation
     * @tparam Vector Output vector type
     * @tparam KeyTransform Key transform type
     * @tparam IndexTransform Index transform type
     * @param result buffer in which to place result
     * @param rhs right hand side buffer
     * @param index_transform transform to apply by index (dense)
     */
    template<typename Vector, typename KeyTransform, typename IndexTransform>
    void
    square_buffered_apply_binary_transform(Vector& result, const dense_vector& rhs, KeyTransform /*key_transform*/,
                                           IndexTransform index_transform) const
    {
        if (empty() || rhs.empty()) {
            return;
        }

        dense_vector& d_result = dtl::vector_base_access::convert(result);
        d_result.resize_to_dimension(std::max(dimension(), rhs.dimension()));

        assert(d_result.m_data.size() >= std::max(dimension(), rhs.dimension()));

        index_transform(&d_result.m_data[0], &m_data[0], &rhs.m_data[0], dimension(), rhs.m_dimension);
    }

public:
    /**
     * @brief  Apply a transform inplace with buffering
     * @tparam Transform Transform type
     * @param result buffer in which to place result (temporarily)
     * @param transform transform to apply
     * @param max_deg Maximum degree
     */
    template<typename Transform>
    void buffered_apply_unary_transform(dense_vector& result, Transform transform, const DEG max_deg) const
    {
        if (empty()) {
            return;
        }

        result.resize_to_dimension(transform.dense_resize(dimension()));
        typename Transform::index_transform it(transform.get_index_transform());

        it(&result.m_data[0], result.dimension(), &m_data[0], dimension(), max_deg);
    }

    /**
     * @brief  Apply a transform inplace with buffering
     * @tparam Transform Transform type
     * @param result buffer in which to place result (temporarily)
     * @param transform transform to apply
     */
    template<typename Transform>
    void buffered_apply_unary_transform(dense_vector& result, Transform transform) const
    {
        if (empty()) {
            return;
        }

        result.resize_to_dimension(transform.dense_resize(dimension()));
        typename Transform::index_transform it(transform.get_index_transform());

        it(&result.m_data[0], result.dimension(), &m_data[0], dimension());
    }

#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
private:
    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive& ar, const unsigned /*version*/)
    {
        ar& m_dimension;
        ar& m_degree;
        ar& m_data;
    }
#endif
};

#undef DECLARE_FUSED_OP

namespace dtl {

template<typename Basis, typename Coeffs>
struct data_access_base<dense_vector<Basis, Coeffs>> {

    using vector_type = dense_vector<Basis, Coeffs>;
    using tag = access_type_dense;

    static const typename Coeffs::S* range_begin(const vector_type& vect)
    {
        return vect.m_data.begin();
    }

    static const typename Coeffs::S* range_end(const vector_type& vect)
    {
        return vect.m_data.end();
    }

    static typename Coeffs::S* range_begin(vector_type& vect)
    {
        return vect.m_data.begin();
    }

    static typename Coeffs::S* range_end(vector_type& vect)
    {
        return vect.m_data.end();
    }



};

}// namespace dtl

}// namespace vectors
}// namespace alg

#endif// LIBALGEBRA_DENSE_VECTOR_H
