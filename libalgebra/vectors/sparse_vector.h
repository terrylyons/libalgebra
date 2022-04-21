/* *************************************************************
Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurk� and Arend Janssen.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */

//  sparse_vector.h

// Include once wrapper
#ifndef DJC_COROPA_LIBALGEBRA_SPARSEVECTORH_SEEN
#define DJC_COROPA_LIBALGEBRA_SPARSEVECTORH_SEEN

#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/unordered_map.hpp>
#endif

#include <iosfwd>
#include <unordered_map>
#include <utility>
#include <vector>

#include "libalgebra/utils/order_trait.h"
#include "libalgebra/vectors/base_vector.h"
#include "libalgebra/vectors/iterators.h"

namespace alg {
namespace vectors {

/// A class to store and manipulate sparse vectors.

// Unordered and Ordered forms
//  sparse_vector is by default ordered (unless the UNORDERED macro is defined)
//  unordered_sparse_vector is not ordered; it is significantly faster for
//  normal functions however iterators may be invalidated by any sort of
//  insertion in the unordered settings
/**
 * @brief Sparse vector where elements are stored as key-value pairs
 *
 * An instance of the sparse_vector class is just a(n unordered) MAP between KEY
 * and SCALAR, with vector space operators. It is a vector of basis elements of
 * type KEY, stored in a MAP class, associated to coefficients given by SCALAR
 * instances. Each basis element refers to the static instance of type BASIS.
 *
 * The MAP class must comes with a std::map<KEY, SCALAR> interface.
 * The scalar type SCALAR corresponds to MAP::mapped_type.
 * By default, the MAP class is taken from the BASIS via the BASIS::MAP
 * type, cf. forward declaration of the sparse_vector class in libalgebra.h.
 *
 * The SCALAR type must come with operators making it an associative
 * algebra (non necessarily commutative) containing the integers (via a
 * suitable constructor). Thus, operators *,+,- must be implemented.
 * It is necessary that the class can be initialized from 0, +1, -1.
 *
 * There is a compatibility condition between the BASIS and MAP classes
 * since the MAP::key_type type and the BASIS::KEY must be the same.
 *
 * For unordered MAP use it is assumed that the follow:
 * References and iterators to the erased elements are invalidated.
 * Other iterators and references are not invalidated. Moreover (C++2014)
 * the internal order of the elements not erased is preserved. However
 * insertion causes a rehash which disrupts all iterators.
 *
 * @tparam Basis Basis of the vector space
 * @tparam Coeffs Coefficient field
 * @tparam MapType The underlying map type in which the data is stored.
 */
template<typename Basis, typename Coeffs, typename MapType = std::unordered_map<typename Basis::KEY, typename Coeffs::S>>
class sparse_vector : /*private*/ MapType, protected base_vector<Basis, Coeffs>
{
    typedef MapType MAP;
    typedef base_vector<Basis, Coeffs> BASE_VEC;

    friend class dtl::data_access_base<sparse_vector>;

public:
    typedef Basis BASIS;
    using MAP::operator[];
    /// Import of set this instance to the zero instance
    using MAP::clear;
    /// Import empty()
    using MAP::empty;
    /// Import size()
    using MAP::size;

    /// Static variables defined in the base vector
    using BASE_VEC::basis;
    using BASE_VEC::mone;
    using BASE_VEC::one;
    using BASE_VEC::zero;

    using BASE_VEC::degree_tag;

    /// Swap the vector instance controlled by *this with the one in the RHS
    void swap(sparse_vector& rhs)
    {
        MAP::swap((MAP&)rhs);
    }

    typedef Coeffs COEFFS;

    /// Import of the KEY type from the MAP class.
    typedef typename BASIS::KEY KEY;
    /// Import of the SCALAR and RATIONAL types from FIELD
    typedef typename COEFFS::S SCALAR;
    typedef typename COEFFS::Q RATIONAL;
    /// Import of the iterator type from the MAP type.
    // typedef typename MAP::iterator iterator;
    /// Import of the KEY constant iterator type from the MAP type.
    // typedef typename MAP::const_iterator const_iterator;

    class iterator_item
    {
        friend class iterators::vector_iterator<iterator_item>;

        friend class sparse_vector;

    public:
        typedef KEY key_type;
        typedef SCALAR& value_type;

        iterator_item()
            : m_iterator()
        {}

        iterator_item(sparse_vector&, typename MAP::iterator it)
            : m_iterator(it)
        {}

        key_type key()
        {
            return m_iterator->first;
        }

        value_type value()
        {
            return m_iterator->second;
        }

        DIMN index() const
        {
            return sparse_vector::basis.key_to_index(m_iterator->first);
        }

        bool operator==(const iterator_item& other) const
        {
            return compare_iterators(other);
        }

        bool operator!=(const iterator_item& other) const
        {
            return !compare_iterators(other);
        }

    private:
        typename MAP::iterator m_iterator;

    private:
        bool compare_iterators(const iterator_item& other) const
        {
            return (m_iterator == other.m_iterator);
        }

        void advance()
        {
            ++m_iterator;
        }
    };

    class const_iterator_item
    {
        friend class iterators::vector_iterator<const_iterator_item>;

        friend class sparse_vector;

    public:
        typedef KEY key_type;
        typedef const SCALAR& value_type;

        const_iterator_item()
            : m_iterator()
        {}

        const_iterator_item(const sparse_vector&, typename MAP::const_iterator it)
            : m_iterator(it)
        {}

        /*
        const_iterator_item& operator=(const const_iterator_item& other)
        {
            m_iterator = other.m_iterator;
            return *this;
        }
    */

        key_type key()
        {
            return m_iterator->first;
        }

        value_type value()
        {
            return m_iterator->second;
        }

        DIMN index() const
        {
            return sparse_vector::basis.key_to_index(m_iterator->first);
        }

        bool operator==(const const_iterator_item& other) const
        {
            return compare_iterators(other);
        }

        bool operator!=(const const_iterator_item& other) const
        {
            return !compare_iterators(other);
        }

    private:
        typename MAP::const_iterator m_iterator;

    private:
        bool compare_iterators(const const_iterator_item& other) const
        {
            return m_iterator == other.m_iterator;
        }

        void advance()
        {
            ++m_iterator;
        }
    };

    typedef iterators::vector_iterator<iterator_item> iterator;
    typedef iterators::vector_iterator<const_iterator_item> const_iterator;

private:
    typename MAP::iterator map_begin()
    {
        return MAP::begin();
    }

    typename MAP::iterator map_end()
    {
        return MAP::end();
    }

    typename MAP::const_iterator map_begin() const
    {
        return MAP::begin();
    }

    typename MAP::const_iterator map_end() const
    {
        return MAP::end();
    }

    typename MAP::iterator map_find(const KEY& key)
    {
        return MAP::find(key);
    }

    typename MAP::const_iterator map_find(const KEY& key) const
    {
        return MAP::find(key);
    }

public:
    // Iterator methods

    /// Iterator to start of vector
    iterator begin()
    {
        if (empty()) {
            return iterator(*this, map_end());
        }
        return iterator(*this, map_begin());
    }

    /// Iterator to end of vector
    iterator end()
    {
        return iterator(*this, map_end());
    }

    /// Const iterator to start of vector
    const_iterator begin() const
    {
        if (empty()) {
            return const_iterator(*this, map_end());
        }
        return const_iterator(*this, map_begin());
    }

    /// Const iterator to end of vector
    const_iterator end() const
    {
        return const_iterator(*this, map_end());
    }

    /// Const iterator to start of vector
    const_iterator cbegin() const
    {
        return begin();
    }

    /// Const iterator to the end of vector
    const_iterator cend() const
    {
        return end();
    }

    /// Get the iterator corresponding to key
    iterator find(const KEY& key)
    {
        return iterator(*this, map_find(key));
    }

    /// Get the const iterator corresponding to key
    const_iterator find(const KEY& key) const
    {
        return const_iterator(*this, map_find(key));
    }

    /// Import the insert from iterator function
    using MAP::insert;

    // Redefine the other inserts
    /// Insert an element from key-value pair
    std::pair<iterator, bool> insert(const std::pair<const KEY, SCALAR>& value)
    {
        if (zero == value.second) {
            return std::pair<iterator, bool>(iterator(*this, MAP::find(value.first)), false);
        }
        std::pair<typename MAP::iterator, bool> p = MAP::insert(value);
        return std::pair<iterator, bool>(iterator(*this, p.first), p.second);
    }

    // Redefine the other inserts
    /// Insert an element from key-value pair
    std::pair<iterator, bool> insert(std::pair<const KEY, SCALAR>& value)
    {
        if (zero == value.second) {
            return std::pair<iterator, bool>(iterator(*this, MAP::find(value.first)), false);
        }
        std::pair<typename MAP::iterator, bool> p = MAP::insert(value);
        return std::pair<iterator, bool>(iterator(*this, p.first), p.second);
    }

    /// Insert an element with using hint to position
    iterator insert(iterator position, const std::pair<const KEY, SCALAR>& value)
    {
        typename MAP::iterator it = MAP::insert(position->m_iterator, value);
        return iterator(*this, it);
    }

    /// Import of erase a KEY from the sparse vector
    using MAP::erase;

    // Redefine the erases involving iterators
    /// Erase the value pointed to by iterator
    void erase(iterator position)
    {
        MAP::erase(position->m_iterator);
    }

    /// Erase a range of elements
    void erase(iterator first, iterator last)
    {
        MAP::erase(first->m_iterator, last->m_iterator);
    }

public:
    /// Given a const instance of a sparse vector, returns a const reference to
    /// the scalar associated to the named basis element. (The default SCALAR
    /// element zero if the basis vector was not present in this sparse vector
    /// instance).
    inline const SCALAR& operator[](const KEY& k) const
    {
        const_iterator found = find(k);
        return (found == cend()) ? zero : found->value();
    }

private:
    void insert_from_pointer(DIMN offset, SCALAR const* begin, SCALAR const* end)
    {
        KEY k(basis.begin());
        while (offset) {
            --offset;
            k = basis.nextkey(k);
        }
        for (SCALAR const* it(begin); it != end && k != basis.end(); ++it, k = basis.nextkey(k)) {
            MAP::insert(std::make_pair(k, *it));
        }
    }

public:
    /*
     * Following the rule of 5, if we define any one of trivial ctor, copy/move ctor or assignment,
     * then we should define all of them. We're going to do this, but let's let the compiler fill
     * them in.
     */

    /** @brief Default constructor.
     * Create an instance of an empty vector.
     * Such a vector is a neutral element for operator+= and operator-=.
     */
    sparse_vector() = default;

    /// Copy constructor.
    sparse_vector(const sparse_vector& v) = default;
    /// Move constructor
    sparse_vector(sparse_vector&& other) noexcept = default;
    /// Copy assignment
    sparse_vector& operator=(const sparse_vector& other) = default;
    /// Move assignment
    sparse_vector& operator=(sparse_vector&& other) noexcept = default;

    /// Unidimensional constructor.
    /**
     * Constructs a sparse_vector corresponding the unique basis
     * element k with coefficient s (+1 by default).
     */
    explicit sparse_vector(const KEY& k, const SCALAR& s = one)
    {
        if (zero != s) {
            (*this)[k] = s;
        }
    }

    /**
     * @brief Construct from pointer to data
     * @param begin Pointer to start of data
     * @param end Pointer to end of data
     */
    sparse_vector(SCALAR const* begin, SCALAR const* end)
        : MAP()
    {
        insert_from_pointer(0, begin, end);
    }

    /**
     * @brief Construct from pointer to data with offset
     *
     * @param offset Offset from beginning of basis order
     * @param pointer to start of data
     * @param end Pointer to end of data
     */
    sparse_vector(DIMN offset, SCALAR const* begin, SCALAR const* end)
        : MAP()
    {
        insert_from_pointer(offset, begin, end);
    }

public:
    /// Returns an instance of the additive inverse of the instance.
    inline sparse_vector operator-() const
    {
        if (empty()) {
            return *this;
        }
        const_iterator in;
        sparse_vector result;
        for (in = begin(); in != end(); ++in) {
            result[in->key()] = -(in->value());
        }
        return result;
    }

    /// Multiplies the instance with scalar s.
    inline sparse_vector& operator*=(const SCALAR& s)
    {
        if (s != zero) {
            iterator it;
            if (!empty()) {
                for (it = begin(); it != end(); ++it) {
                    it->value() *= s;
                }
            }
        }
        else {
            clear();
        }
        return *this;
    }

    /// Binary version of operator*=()
    inline sparse_vector operator*(const SCALAR& rhs) const
    {
        sparse_vector result(*this);
        return result *= rhs;
    };

    /// Divides the instance by scalar s.
    inline sparse_vector& operator/=(const RATIONAL& s)
    {
        iterator it;
        if (!empty()) {
            for (it = begin(); it != end(); ++it) {
                RATIONAL temp(1);
                it->value() *= (temp / s);
            }
        }
        return *this;
    }

    /// Binary instance of  operator/=()
    inline sparse_vector operator/(const RATIONAL& rhs) const
    {
        sparse_vector result(*this);
        return result /= rhs;
    };

    /// Adds a sparse_vector to the instance.
    inline sparse_vector& operator+=(const sparse_vector& rhs)
    {
        iterator it;
        const_iterator cit;
        if (rhs.empty()) {
            return *this;
        }
        if (empty()) {
            return *this = rhs;
        }
        for (cit = rhs.begin(); cit != rhs.end(); ++cit) {// Instead of a bare (*this)[cit->first] += cit->second;
            it = find(cit->key());
            if (it == end()) {
                (*this)[cit->key()] = cit->value();
            }
            else if ((it->value() += cit->value()) == zero) {
                erase(it->key());
            }
        }
        return *this;
    }

    /// Binary version of  operator+=()
    inline sparse_vector operator+(const sparse_vector& rhs) const
    {
        sparse_vector result(*this);
        result += rhs;
        return result;
    };

    /// Subtracts a sparse_vector to the instance.
    inline sparse_vector& operator-=(const sparse_vector& rhs)
    {
        iterator it;
        const_iterator cit;
        if (rhs.empty()) {
            return *this;
        }
        if (empty()) {
            return *this = -rhs;
        }
        for (cit = rhs.begin(); cit != rhs.end(); ++cit) {// Instead of a bare (*this)[cit->first] -= cit->second;
            it = find(cit->key());
            if (it == end()) {
                (*this)[cit->key()] = -(cit->value());
            }
            else if ((it->value() -= cit->value()) == zero) {
                erase(it->key());
            }
        }
        return *this;
    }

    /// Binary version of  operator-=()
    inline sparse_vector operator-(const sparse_vector& rhs) const
    {
        sparse_vector result(*this);
        result -= rhs;
        return result;
    };

private:
    void min_impl(const sparse_vector& rhs, std::true_type)
    {
        typename MAP::iterator it(map_begin()), itend(map_end());
        typename MAP::const_iterator cit(rhs.map_begin()), cend(rhs.map_end());
        for (; it != itend && cit != cend;) {
            int c = (it->first < cit->first)    ? 1
                    : (cit->first < it->first)  ? 2
                    : (cit->first == it->first) ? 3
                                                : 4;
            switch (c) {
            case 1: {
                if (it->second >= zero)
                    erase(it++);
                break;
            }
            case 2: {
                if (cit->second < zero)
                    insert(*cit);
                ++cit;
                break;
            }
            case 3: {
                operator[](it->first) =
                        ((it->second < cit->second) ? (it->second) : (cit->second));
                ++cit;
                ++it;
                break;
            }
            default:;
            }
        }
        if (cit == cend) {
            for (; it != itend;)
                if (it->second >= zero)
                    erase(it++);
                else
                    ++it;
        }
        if (it == itend) {
            for (; cit != cend; ++cit)
                if (cit->second < zero)
                    insert(*cit);
        }
    }

    void min_impl(const sparse_vector& rhs, std::false_type)
    {

        typename std::vector<std::pair<KEY, SCALAR>> target(map_begin(), map_end()), source(rhs.map_begin(), rhs.map_end());
        const auto& comp = [](typename std::pair<KEY, SCALAR> lhs, typename std::pair<KEY, SCALAR> rhs) -> bool {
            return lhs.first < rhs.first;
        };
        std::sort(target.begin(), target.end(), comp);
        std::sort(source.begin(), source.end(), comp);
        typename std::vector<std::pair<KEY, SCALAR>>::iterator it = target.begin();
        typename std::vector<std::pair<KEY, SCALAR>>::const_iterator cit = source.begin();
        for (; it != target.end() && cit != source.end();) {
            int c = (it->first < cit->first) ? 1 : (cit->first < it->first) ? 2
                    : (cit->first == it->first)                             ? 3
                                                                            : 4;
            switch (c) {
            case 1: {
                if (it->second >= zero) {
                    erase((it++)->first);
                }
                break;
            }
            case 2: {
                if (cit->second < zero) {
                    insert(*cit);
                }
                ++cit;
                break;
            }
            case 3: {
                operator[](it->first) = ((it->second < cit->second) ? (it->second) : (cit->second));
                ++cit;
                ++it;
                break;
            }
            default:;
            }
        }
        if (cit == source.end()) {
            for (; it != target.end();) {
                if (it->second >= zero) {
                    erase((it++)->first);
                }
                else {
                    ++it;
                }
            }
        }
        if (it == target.end()) {
            for (; cit != source.end(); ++cit) {
                if (cit->second < zero) {
                    insert(*cit);
                }
            }
        }
    }

public:
    /// Where SCA admits an order forms the min of two sparse vectors
    inline sparse_vector& operator&=(const sparse_vector& rhs)
    {
        // these min max operators are slower (factor of 3?) on unordered sparse vectors
        min_impl(rhs, utils::is_ordered<MAP>());
        return *this;
    }

    /// Binary version of  operator&=()
    inline sparse_vector operator&(const sparse_vector& rhs) const
    {
        sparse_vector result(*this);
        result &= rhs;
        return result;
    };

private:
    void max_impl(const sparse_vector& rhs, std::true_type)
    {
        typename MAP::iterator it(map_begin()), itend(map_end());
        typename MAP::const_iterator cit(rhs.map_begin()), cend(rhs.map_end());
        for (; it != itend && cit != cend;) {
            // c++11 syntax auto
            int c = (it->first < cit->first)    ? 1
                    : (cit->first < it->first)  ? 2
                    : (cit->first == it->first) ? 3
                                                : 4;
            switch (c) {
            case 1: {
                if (it->second <= zero)
                    erase(it++);
                break;
            }
            case 2: {
                if (cit->second > zero)
                    insert(*cit);
                ++cit;
                break;
            }
            case 3: {
                operator[](it->first) =
                        ((it->second > cit->second) ? (it->second) : (cit->second));
                ++cit;
                ++it;
                break;
            }
            default:;
            }
        }
        if (cit == cend) {
            for (; it != itend;)
                if (it->second <= zero)
                    erase(it++);
                else
                    ++it;
        }
        if (it == itend) {
            for (; cit != cend; ++cit)
                if (cit->second > zero)
                    insert(*cit);
        }
    }

    /// Where SCA admits an order forms the max of two sparse vectors
    void max_impl(const sparse_vector& rhs, std::false_type)
    {
        typename std::vector<std::pair<KEY, SCALAR>> target(map_begin(), map_end()), source(rhs.map_begin(), rhs.map_end());
        std::sort(target.begin(), target.end(), comp);
        std::sort(source.begin(), source.end(), comp);

        typename std::vector<std::pair<KEY, SCALAR>>::iterator it = target.begin();
        typename std::vector<std::pair<KEY, SCALAR>>::const_iterator cit = source.begin();
        for (; it != target.end() && cit != source.end();) {
            auto c = (it->first < cit->first) ? 1 : (cit->first < it->first) ? 2
                    : (cit->first == it->first)                              ? 3
                                                                             : 4;
            switch (c) {
            case 1: {
                if (it->second <= zero) {
                    erase((it++)->first);
                }
                break;
            }
            case 2: {
                if (cit->second > zero) {
                    insert(*cit);
                }
                ++cit;
                break;
            }
            case 3: {
                operator[](it->first) = ((it->second > cit->second) ? (it->second) : (cit->second));
                ++cit;
                ++it;
                break;
            }
            default:;
            }
        }
        if (cit == source.end()) {
            for (; it != target.end();) {
                if (it->second <= zero) {
                    erase((it++)->first);
                }
                else {
                    ++it;
                }
            }
        }
        if (it == target.end()) {
            for (; cit != source.end(); ++cit) {
                if (cit->second > zero) {
                    insert(*cit);
                }
            }
        }
    }

public:
    /// Where SCA admits an order forms the max of two sparse vectors
    inline sparse_vector& operator|=(const sparse_vector& rhs)
    {
        max_impl(rhs, utils::is_ordered<MAP>());
        return *this;
    }

    /// Binary version of  operator|=()
    inline sparse_vector operator|(const sparse_vector& rhs) const
    {
        sparse_vector result(*this);
        result |= rhs;
        return result;
    };

    /// A version of operator+=(rhs.scal_prod(s))
    /// when RHS is a scaled basis vector
    inline sparse_vector& add_scal_prod(const KEY& rhs, const SCALAR& s)
    {
        // sparse addition
        if (zero == (operator[](rhs) += s)) {
            erase(rhs);
        }
        return *this;
    }

    /// A version of operator+=(rhs.scal_prod(s))
    /// when RHS is a sparse vector scaled on right
    inline sparse_vector& add_scal_prod(const sparse_vector& rhs, const SCALAR& s)
    {
        if ((s == zero) || rhs.empty()) {
            return *this;
        }
        if (empty()) {
            *this = rhs;
            return operator*=(s);
        }
        iterator it = begin();
        const_iterator cit;
        for (cit = rhs.begin(); cit != rhs.end(); ++cit) {          // Instead of a bare (*this)[cit->first] += cit->second * s;
            it = this->insert(it, std::make_pair(cit->key(), zero));// note this fails if the entry is already
            // there but sets it in any case
            if ((it->value() += cit->value() * s) == zero)
            // erase returns void until c++11
            {
                iterator j(it++);
                erase(j);
            }
            else {
                ++it;
            }
        }
        return *this;
    }

    /// A version of operator-=(rhs.scal_prod(s))
    /// when RHS is a scaled basis vector
    inline sparse_vector& sub_scal_prod(const KEY& rhs, const SCALAR& s)
    {
        // sparse addition
        if (zero == (operator[](rhs) -= s)) {
            erase(rhs);
        }
        return *this;
    }

    /// A version of operator-=(rhs.scal_prod(s))
    /// when RHS is a sparse vector scaled on right
    inline sparse_vector& sub_scal_prod(const sparse_vector& rhs, const SCALAR& s)
    {
        iterator it;
        const_iterator cit;
        if ((s == zero) || rhs.empty()) {
            return *this;
        }
        if (empty()) {
            *this = rhs;
            return operator*=(-s);
        }
        for (cit = rhs.begin(); cit != rhs.end(); ++cit) {// Instead of a bare (*this)[cit->first] -= cit->second * s;
            it = find(cit->key());
            if (it == end()) {
                (*this)[cit->key()] = cit->value() * -s;
            }
            else if ((it->value() -= cit->value() * s) == zero) {
                erase(it->key());
            }
        }
        return *this;
    }

    /// A fast version of operator+=(rhs.scal_div(s))
    inline sparse_vector& add_scal_div(const sparse_vector& rhs, const RATIONAL& s)
    {
        iterator it;
        const_iterator cit;
        if (rhs.empty()) {
            return *this;
        }
        if (empty()) {
            *this = rhs;
            return operator/=(s);
        }
        for (cit = rhs.begin(); cit != rhs.end(); ++cit) {// Instead of a bare (*this)[cit->first] += cit->second / s;
            it = find(cit->key());
            if (it == end()) {
                (*this)[cit->key()] = cit->value() / s;
            }
            else if ((it->value() += (cit->value() / s)) == zero) {
                erase(it->key());
            }
        }
        return *this;
    }

    /// A fast version of operator-=(rhs.scal_div(s))
    inline sparse_vector& sub_scal_div(const sparse_vector& rhs, const RATIONAL& s)
    {
        iterator it;
        const_iterator cit;
        if (rhs.empty()) {
            return *this;
        }
        if (empty()) {
            *this = rhs;
            return operator/=(-s);
        }
        for (cit = rhs.begin(); cit != rhs.end(); ++cit) {// Instead of a bare (*this)[cit->first] -= cit->second / s;
            it = find(cit->key());
            if (it == end()) {
                (*this)[cit->key()] = -cit->value() / s;
            }
            else if ((it->value() -= (cit->value() / s)) == zero) {
                erase(it->key());
            }
        }
        return *this;
    }

    /// Fused inplace addition with scalar division
    inline sparse_vector& add_scal_div(const KEY& rhs, const RATIONAL& s)
    {
        if (zero == (operator[](rhs) += one / s)) {
            erase(rhs);
        }
        return *this;
    }

    /// Fused inplace subtraction with rational division
    inline sparse_vector& sub_scal_div(const KEY& rhs, const RATIONAL& s)
    {
        if (zero == (operator[](rhs) -= one / s)) {
            erase(rhs);
        }
        return *this;
    }

    /// Compares the instance to a sparse_vector.
    bool operator==(const sparse_vector& rhs) const
    {
        if (size() != rhs.size()) {
            return false;
        }
        const_iterator i, j, jend(rhs.end()), iend(end());
        for (i = begin(); i != iend; ++i) {
            j = rhs.find(i->key());
            if ((j == jend) || (j->value() != i->value())) {
                return false;
            }
        }
        return true;
    }

    /// Lexicographically compares the instance to a sparse_vector.
    bool operator<(const sparse_vector& rhs) const
    {
        return std::lexicographical_compare(map_begin(), map_end(), rhs.map_begin(), rhs.map_end());
    }

    /// Boolean negation of operator==()
    bool operator!=(const sparse_vector& rhs) const
    {
        return !operator==(rhs);
    }

    /// Get the current maximum degree of elements in the vector
    DEG degree() const
    {
        DEG ans(0);
        for (const_iterator it(begin()); it != end(); ++it) {
            ans = std::max(basis.degree(it->key()), ans);
        }
        return ans;
    }

    /// Check if the current maximum degree is equal to given value
    bool degree_equals(const DEG degree) const
    {
        bool result(false);
        DEG d;
        for (const_iterator it(begin()); it != end(); ++it) {
            d = basis.degree(it->key());
            if (d > degree) {
                return false;
            }
            else if (!result && d == degree) {
                result = true;
            }
        }
        return result;
    }

    /// Computes the l1 norm of sparse vector with respect to this basis
    inline SCALAR NormL1() const
    {
        const_iterator i;
        SCALAR ans(zero);
        for (i = begin(); i != end(); ++i) {
            ans += abs(i->value());
        }
        return ans;
    }

    /// Computes the l1 norm of degree d component of a sparse vector with respect
    /// to this basis
    inline SCALAR NormL1(const DEG& d) const
    {
        const_iterator i;
        SCALAR ans(zero);
        for (i = begin(); i != end(); ++i) {
            if (d == basis.degree(i->key())) {
                ans += abs(i->value());
            }
        }
        return ans;
    }

    /// Compute the l-infinity norm of a sparse vector
    SCALAR NormLInf() const
    {
        const_iterator i;
        SCALAR ans(zero);
        for (i = begin(); i != end(); ++i) {
            ans = std::max(abs(i->value()), ans);
        }
        return ans;
    }

    /// Computes the l-infinity norm of degree d component of a sparse vector with
    /// respect to this basis
    inline SCALAR NormLInf(const DEG& d) const
    {
        const_iterator i;
        SCALAR ans(zero);
        for (i = begin(); i != end(); ++i) {
            if (d == basis.degree(i->key())) {
                ans = std::max(abs(i->value()), ans);
            }
        }
        return ans;
    }

    /// Outputs a sparse_vector to an std::ostream.
    /**
    It is assumed that there is std::ostream support to
    output std::pair<BASIS*, KEY> and SCALAR types.
    */

    static bool comp(typename std::pair<KEY, SCALAR> lhs, typename std::pair<KEY, SCALAR> rhs)
    {
        return lhs.first < rhs.first;
    };

private:
    void print_members_impl(std::ostream& os, std::false_type) const
    {
        std::pair<BASIS*, KEY> token;
        token.first = &sparse_vector::basis;
        typename std::vector<std::pair<KEY, SCALAR>>::const_iterator cit;
        typename std::vector<std::pair<KEY, SCALAR>> buffer(map_begin(), map_end());
        std::sort(buffer.begin(), buffer.end(), comp);
        for (cit = buffer.begin(); cit != buffer.end(); ++cit) {
            token.second = cit->first;
            os << ' ' << cit->second << '(' << token << ')';
        }
    }

    void print_members_impl(std::ostream& os, std::true_type) const
    {
        std::pair<BASIS*, KEY> token;
        token.first = &sparse_vector::basis;
        const_iterator cit;

        for (cit = begin(); cit != end(); ++cit) {
            token.second = cit->key();
            os << ' ' << cit->value() << '(' << token << ')';
        }
    }

protected:
    /// Print members to output stream
    void print_members(std::ostream& os) const
    {
        print_members_impl(os, utils::is_ordered<MAP>());
    }

public:
    /// Stream out operator
    inline friend std::ostream& operator<<(std::ostream& os, const sparse_vector& rhs)
    {
        os << '{';
        rhs.print_members(os);
        os << " }";
        return os;
    }

protected:
    void fill_buffer(std::vector<std::pair<KEY, SCALAR>>& buffer, std::true_type) const
    {
        buffer.assign(map_begin(), map_end());
    }

    void fill_buffer(std::vector<std::pair<KEY, SCALAR>>& buffer, std::false_type) const
    {
        buffer.assign(map_begin(), map_end());
        std::sort(buffer.begin(), buffer.end(), typename BASIS::ordering_tag::pair_order());
    }

    /// copy the (key, value) elements from rhs to a sorted vector buffer (using
    /// the key for sorting) and construct an increasing vector iterators so that
    /// segment [iterators[i-1], iterators[i]) contains keys of degree i; the
    /// first begins at [begin(), and the last ends at end), and it can be empty
    void separate_by_degree(std::vector<std::pair<KEY, SCALAR>>& buffer, const sparse_vector& rhs, const size_t DEPTH1,
                            std::vector<typename std::vector<std::pair<KEY, SCALAR>>::const_iterator>& iterators) const
    {
        rhs.fill_buffer(buffer, utils::is_ordered<MAP>());

        iterators.assign(DEPTH1 + 1, buffer.end());
        unsigned deg = 0;
        for (typename std::vector<std::pair<KEY, SCALAR>>::const_iterator j0 = buffer.begin(); j0 != buffer.end();
             j0++) {
            DEG d = basis.degree(j0->first);
            assert(d >= deg && d <= DEPTH1);// order assumed to respect degree
            while (deg < d) {
                iterators[deg++] = j0;
            }
            // deg == d
        }
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
    triangular_buffered_apply_binary_transform(Vector& result, const sparse_vector& rhs, KeyTransform key_transform,
                                               const DEG /* max_depth*/) const
    {
        // Unused parameter max_depth?
        //DEG max_degree = std::min(max_depth, degree_tag.max_degree);
        DEG max_degree = degree_tag.max_degree;

        // create buffers to avoid unnecessary calls to MAP inside loop
        std::vector<std::pair<KEY, SCALAR>> buffer;
        std::vector<typename std::vector<std::pair<KEY, SCALAR>>::const_iterator> iterators;
        separate_by_degree(buffer, rhs, max_degree, iterators);

        typename std::vector<std::pair<KEY, SCALAR>>::const_iterator j;
        const_iterator i(begin()), iEnd(end());
        DEG rhdegree;
        for (; i != iEnd; ++i) {
            const KEY& k = i->key();
            rhdegree = max_degree - basis.degree(k);
            typename std::vector<std::pair<KEY, SCALAR>>::const_iterator& jEnd = iterators[rhdegree];
            for (j = buffer.begin(); j != jEnd; ++j) {
                key_transform(result, k, i->value(), j->first, j->second);
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
    triangular_buffered_apply_binary_transform(Vector& result, const sparse_vector& rhs, KeyTransform key_transform,
                                               IndexTransform /* index_transform */, const DEG max_depth) const
    {
        triangular_buffered_apply_binary_transform(result, rhs, key_transform, max_depth);
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
    square_buffered_apply_binary_transform(Vector& result, const sparse_vector& rhs, KeyTransform key_transform) const
    {
        // create buffer to avoid unnecessary calls to MAP inside loop
        std::vector<std::pair<KEY, SCALAR>> buffer(rhs.map_begin(), rhs.map_end());
        const_iterator i;

        // DEPTH1 == 0
        typename std::vector<std::pair<KEY, SCALAR>>::const_iterator j;
        for (i = begin(); i != end(); ++i) {
            for (j = buffer.begin(); j != buffer.end(); ++j) {
                key_transform(result, i->key(), i->value(), j->first, j->second);
            }
        }
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
    triangular_unbuffered_apply_binary_transform(const sparse_vector& rhs, KeyTransform key_transform,
                                                 IndexTransform index_transform, const DEG max_depth)
    {
        sparse_vector result;
        triangular_buffered_apply_binary_transform(result, rhs, key_transform, index_transform, max_depth);
        swap(result);
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
    triangular_unbuffered_apply_binary_transform(const sparse_vector& rhs, KeyTransform key_transform,
                                                 const DEG max_depth)
    {
        sparse_vector result;
        triangular_buffered_apply_binary_transform(result, rhs, key_transform, max_depth);
        swap(result);
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
    square_buffered_apply_binary_transform(Vector& result, const sparse_vector& rhs, KeyTransform key_transform,
                                           IndexTransform /* index_transform */) const
    {
        square_buffered_apply_binary_transform(result, rhs, key_transform);
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
    void buffered_apply_unary_transform(sparse_vector& result, Transform transform, const DEG max_deg) const
    {
        if (empty()) {
            return;
        }

        typename Transform::key_transform kt(transform.get_key_transform());
        kt(result, *this, max_deg);
    }

    /**
     * @brief  Apply a transform inplace with buffering
     * @tparam Transform Transform type
     * @param result buffer in which to place result (temporarily)
     * @param transform transform to apply
     */
    template<typename Transform>
    void buffered_apply_unary_transform(sparse_vector& result, Transform transform) const
    {
        if (empty()) {
            return;
        }

        typename Transform::key_transform kt(transform.get_key_transform());
        kt(result, *this);
    }
#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
private:
    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive& ar, unsigned int const /* version */)
    {
        ar& boost::serialization::base_object<MapType>(*this);
    }
#endif
};

namespace dtl {

template<typename Basis, typename Coeffs, typename Map>
struct data_access_base<sparse_vector<Basis, Coeffs, Map>> {

    using vector_type = sparse_vector<Basis, Coeffs, Map>;
    using tag = access_type_sparse;

    static typename Map::const_iterator range_begin(const vector_type& vect)
    {
        return vect.map_begin();
    }

    static typename Map::const_iterator range_end(const vector_type& vect)
    {
        return vect.map_end();
    }

    static typename Map::iterator range_begin(vector_type& vect)
    {
        return vect.map_begin();
    }

    static typename Map::iterator range_end(vector_type& vect)
    {
        return vect.map_end();
    }


};

}// namespace dtl

}// namespace vectors
}// namespace alg

// Include once wrapper
// DJC_COROPA_LIBALGEBRA_SPARSEVECTORH_SEEN
#endif
// EOF.
