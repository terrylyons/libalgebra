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

#include "libalgebra/vectors/base_vector.h"
#include "libalgebra/vectors/iterators.h"

namespace alg {
namespace vectors {


// This is a macro because template aliases are c++11
#define LIBALGEBRA_DEFAULT_MAP_TYPE \
    std::map<typename Basis::KEY, typename Coeffs::S>

/// A class to store and manipulate sparse vectors.

//Unordered and Ordered forms
// sparse_vector is by default ordered (unless the UNORDERED macro is defined)
// unordered_sparse_vector is not ordered; it is significantly faster for normal functions
// however iterators may be invalidated by any sort of insertion in the unordered settings
/**
An instance of the sparse_vector class is just a(n unordered) MAP between KEY and
SCALAR, with vector space operators. It is a vector of basis elements
of type KEY, stored in a MAP class, associated to coefficients given by
SCALAR instances. Each basis element refers to the static instance of type
BASIS.

The MAP class must comes with a std::map<KEY, SCALAR> interface.
The scalar type SCALAR correponds to MAP::mapped_type.
By default, the MAP class is taken from the BASIS via the BASIS::MAP
type, cf. forward declaration of the sparse_vector class in libalgebra.h.

The SCALAR type must come with operators making it an associative
algebra (non necessarily commutative) containing the integers (via a
suitable constructor). Thus, operators *,+,- must be implemented.
It is necessary that the class can be initialized from 0, +1, -1.

There is a compatibility condition between the BASIS and MAP classes
since the MAP::key_type type and the BASIS::KEY must be the same.

For unordered MAP use it is assumed that the follow:
References and iterators to the erased elements are invalidated.
Other iterators and references are not invalidated. Moreover (C++2014)
the internal order of the elements not erased is preserved. However
insertion causes a rehash which disrupts all iterators
*/
template<typename Basis, typename Coeffs,
        typename MapType=LIBALGEBRA_DEFAULT_MAP_TYPE >
class sparse_vector : /*private*/ MapType, protected base_vector<Basis, Coeffs>
{
    typedef MapType MAP;
    typedef Basis BASIS;
    typedef base_vector<Basis, Coeffs> BASE_VEC;
public:

    using MAP::operator[];
    /// Import of set this instance to the zero instance
    using MAP::clear;
    /// Import empty()
    using MAP::empty;
    /// Import size()
    using MAP::size;

    /// Static variables defined in the base vector
    using BASE_VEC::basis;
    using BASE_VEC::one;
    using BASE_VEC::mone;
    using BASE_VEC::zero;

    using BASE_VEC::degree_tag;

    /// Swap the vector instance controlled by *this with the one in the RHS
    void swap(sparse_vector &rhs)
    {
        MAP::swap((MAP &) rhs);
    }

    typedef Coeffs FIELD;

    /// Import of the KEY type from the MAP class.
    typedef typename BASIS::KEY KEY;
    /// Import of the SCALAR and RATIONAL types from FIELD
    typedef typename FIELD::S SCALAR;
    typedef typename FIELD::Q RATIONAL;
    /// Import of the iterator type from the MAP type.
    //typedef typename MAP::iterator iterator;
    /// Import of the KEY constant iterator type from the MAP type.
    //typedef typename MAP::const_iterator const_iterator;

    class iterator_item
    {
        friend class iterators::vector_iterator<iterator_item>;

        friend class sparse_vector;

    public:

        typedef KEY key_type;
        typedef SCALAR &value_type;

        iterator_item() : m_iterator()
        {}

        iterator_item(const iterator_item &other)
                : m_iterator(other.m_iterator)
        {}

        iterator_item(sparse_vector &, typename MAP::iterator it)
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

        bool operator==(const iterator_item &other) const
        {
            return compare_iterators(other);
        }

        bool operator!=(const iterator_item &other) const
        {
            return !compare_iterators(other);
        }

    private:
        typename MAP::iterator m_iterator;

    private:

        bool compare_iterators(const iterator_item &other) const
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
        typedef const SCALAR &value_type;

        const_iterator_item() : m_iterator()
        {}

        const_iterator_item(const const_iterator_item &other)
                : m_iterator(other.m_iterator)
        {}

        const_iterator_item(const sparse_vector &,
                            typename MAP::const_iterator it)
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

        bool operator==(const const_iterator_item &other) const
        {
            return compare_iterators(other);
        }

        bool operator!=(const const_iterator_item &other) const
        {
            return !compare_iterators(other);
        }

    private:
        typename MAP::const_iterator m_iterator;

    private:

        bool compare_iterators(const const_iterator_item &other) const
        {
            bool result = m_iterator != other.m_iterator;
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

    typename MAP::iterator map_find(const KEY &key)
    {
        return MAP::find(key);
    }

    typename MAP::const_iterator map_find(const KEY &key) const
    {
        return MAP::find(key);
    }

public:

    // Iterator methods

    iterator begin()
    {
        if (empty()) {
            return iterator(*this, map_end());
        }
        return iterator(*this, map_begin());
    }

    iterator end()
    {
        return iterator(*this, map_end());
    }

    const_iterator begin() const
    {
        if (empty()) {
            return const_iterator(*this, map_end());
        }
        return const_iterator(*this, map_begin());
    }

    const_iterator end() const
    {
        return const_iterator(*this, map_end());
    }

    const_iterator cbegin() const
    {
        return begin();
    }

    const_iterator cend() const
    { return end(); }

    iterator find(const KEY &key)
    {
        return iterator(*this, map_find(key));
    }

    const_iterator find(const KEY &key) const
    {
        return const_iterator(*this, map_find(key));
    }

    /// Import the insert from iterator function
    using MAP::insert;

    // Redefine the other inserts
    std::pair<iterator, bool> insert(const std::pair<const KEY, SCALAR> &value)
    {
        if (zero == value.second) {
            return std::pair<iterator, bool>(
                    iterator(*this, MAP::find(value.first)), false);
        }
        std::pair<typename MAP::iterator, bool> p = MAP::insert(value);
        return std::pair<iterator, bool>(iterator(*this, p.first), p.second);
    }

    // Redefine the other inserts
    std::pair<iterator, bool> insert(std::pair<const KEY, SCALAR> &value)
    {
        if (zero == value.second) {
            return std::pair<iterator, bool>(
                    iterator(*this, MAP::find(value.first)), false);
        }
        std::pair<typename MAP::iterator, bool> p = MAP::insert(value);
        return std::pair<iterator, bool>(iterator(*this, p.first), p.second);
    }

    iterator insert(iterator position, const std::pair<const KEY, SCALAR> &value)
    {
        typename MAP::iterator it = MAP::insert(position->m_iterator, value);
        return iterator(*this, it);
    }

    /// Import of erase a KEY from the sparse vector
    using MAP::erase;

    // Redefine the erases involving iterators
    void erase(iterator position)
    {
        MAP::erase(position->m_iterator);
    }

    void erase(iterator first, iterator last)
    {
        MAP::erase(
                first->m_iterator,
                last->m_iterator
        );
    }


public:

    /// Given a const instance of a sparse vector, returns a const reference to the scalar associated to the named basis element. (The default SCALAR element zero if the basis vector was not present in this sparse vector instance).
    inline const SCALAR &operator[](const KEY& k) const
    {
        const_iterator found = find(k);
        return (found == cend()) ? zero : found->value();
    }


public:
    /// Default constructor.
    /**
    * Create an instance of an empty vector.
    * Such a vector is a neutral element for operator+= and operator-=.
    */
    sparse_vector(void)
    {}

    /// Copy constructor.
    sparse_vector(const sparse_vector &v) : MAP((const MAP &) v)
    {}

    /// Unidimensional constructor.
    /**
    * Constructs a sparse_vector corresponding the unique basis
    * element k with coefficient s (+1 by default).
    */
    explicit sparse_vector(const KEY &k, const SCALAR &s = one)
    {
        if (zero != s)
            (*this)[k] = s;
    }


public:


    /// Returns an instance of the additive inverse of the instance.
    inline sparse_vector operator-(void) const
    {
        if (empty())
            return *this;
        const_iterator in;
        sparse_vector result;
        for (in = begin(); in != end(); ++in)
            result[in->key()] = -(in->value());
        return result;
    }

    /// Multiplies the instance with scalar s.
    inline sparse_vector &operator*=(const SCALAR &s)
    {
        if (s != zero) {
            iterator it;
            if (!empty())
                for (it = begin(); it != end(); ++it)
                    it->value() *= s;
        } else
            clear();
        return *this;
    }

    /// Binary version of operator*=()
    inline __DECLARE_BINARY_OPERATOR(sparse_vector, *, *=, SCALAR);

    /// Divides the instance by scalar s.
    inline sparse_vector &operator/=(const RATIONAL &s)
    {
        iterator it;
        if (!empty())
            for (it = begin(); it != end(); ++it) {
                RATIONAL temp(1);
                it->value() *= (temp / s);
            }
        return *this;
    }

    /// Binary instance of  operator/=()
    inline __DECLARE_BINARY_OPERATOR(sparse_vector,
                                     /, /=, RATIONAL);

    /// Adds a sparse_vector to the instance.
    inline sparse_vector &operator+=(const sparse_vector &rhs)
    {
        iterator it;
        const_iterator cit;
        if (rhs.empty())
            return *this;
        if (empty())
            return *this = rhs;
        for (cit = rhs.begin(); cit != rhs.end(); ++cit) { // Instead of a bare (*this)[cit->first] += cit->second;
            it = find(cit->key());
            if (it == end())
                (*this)[cit->key()] = cit->value();
            else if ((it->value() += cit->value()) == zero)
                erase(it->key());
        }
        return *this;
    }

    /// Binary version of  operator+=()
    inline __DECLARE_BINARY_OPERATOR(sparse_vector,
                                     +, +=, sparse_vector);

    /// Subtracts a sparse_vector to the instance.
    inline sparse_vector &operator-=(const sparse_vector &rhs)
    {
        iterator it;
        const_iterator cit;
        if (rhs.empty())
            return *this;
        if (empty())
            return *this = -rhs;
        for (cit = rhs.begin(); cit != rhs.end(); ++cit) { // Instead of a bare (*this)[cit->first] -= cit->second;
            it = find(cit->key());
            if (it == end())
                (*this)[cit->key()] = -(cit->value());
            else if ((it->value() -= cit->value()) == zero)
                erase(it->key());
        }
        return *this;
    }

    /// Binary version of  operator-=()
    inline __DECLARE_BINARY_OPERATOR(sparse_vector,
                                     -, -=, sparse_vector);

    /// Where SCA admits an order forms the min of two sparse vectors
    inline sparse_vector &operator&=(const sparse_vector &rhs)
    {
// these min max operators are slower (factor of 3?) on unordered sparse vectors
#ifdef UNORDEREDMAP
        {
        typename std::vector<std::pair<KEY, SCALAR> >
                target(map_begin(), map_end()),
                source(rhs.map_begin(), rhs.map_end());
        const auto & comp = [](typename std::pair<KEY, SCALAR>  lhs, typename std::pair<KEY, SCALAR>  rhs)->bool {return lhs.first < rhs.first; };
        std::sort(target.begin(), target.end(), comp);
        std::sort(source.begin(), source.end(), comp);
        typename std::vector<std::pair<KEY, SCALAR> >::iterator it = target.begin();
        typename std::vector<std::pair<KEY, SCALAR> >::const_iterator cit = source.begin();
        for (; it != target.end() && cit != source.end(); )
        {
            int c = (it->first < cit->first) ? 1 : (cit->first < it->first) ? 2 : (cit->first == it->first) ? 3 : 4;
            switch (c)
            {
            case 1: {
                if (!(it->second < SCALAR(0))) erase((it++)->first);
                break;
            }
            case 2: {
                if (cit->second < SCALAR(0)) insert(*cit);
                ++cit;
                break;
            }
            case 3: {
                operator[](it->first) = ((it->second < cit->second) ? (it->second) : (cit->second));
                ++cit;
                ++it;
                break; }
            default:;
            }
        }
        if (cit == source.end())
        {
            for (; it != target.end();)
                if (!(it->second < SCALAR(0))) erase((it++)->first);
                else ++it;
        }
        if (it == target.end())
        {
            for (; cit != source.end(); ++cit)
                if (cit->second < SCALAR(0)) insert(*cit);
        }
    }
#else
        typename MAP::iterator it(map_begin()), itend(map_end());
        typename MAP::const_iterator cit(rhs.map_begin()), cend(rhs.map_end());
        for (; it != itend && cit != cend;) {
            int c = (it->first < cit->first) ? 1 : (cit->first < it->first) ? 2 : (cit->first == it->first) ? 3 : 4;
            switch (c) {
                case 1: {
                    if (!(it->second < SCALAR(0))) erase(it++);
                    break;
                }
                case 2: {
                    if (cit->second < SCALAR(0)) insert(*cit);
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
        if (cit == cend) {
            for (; it != itend;)
                if (!(it->second < SCALAR(0))) erase(it++);
                else ++it;
        }
        if (it == itend) {
            for (; cit != cend; ++cit)
                if (cit->second < SCALAR(0)) insert(*cit);
        }
#endif
        return *this;
    }

    /// Binary version of  operator&=()
    inline __DECLARE_BINARY_OPERATOR(sparse_vector, &, &=, sparse_vector);

    /// Where SCA admits an order forms the max of two sparse vectors
    inline sparse_vector &operator|=(const sparse_vector &rhs)
    {
#ifdef UNORDEREDMAP

        typename std::vector<std::pair<KEY, SCALAR> >
                target(map_begin(), map_end()),
                source(rhs.map_begin(), rhs.map_end());
    std::sort(target.begin(), target.end(), comp);
    std::sort(source.begin(), source.end(), comp);

    typename std::vector<std::pair<KEY, SCALAR> >::iterator it = target.begin();
    typename std::vector<std::pair<KEY, SCALAR> >::const_iterator cit = source.begin();
    for (; it != target.end() && cit != source.end(); )

    {
        auto c = (it->first < cit->first) ? 1 : (cit->first < it->first) ? 2 : (cit->first == it->first) ? 3 : 4;
        switch (c)
        {
        case 1: {
            if (!(it->second > SCALAR(0))) erase((it++)->first);
            break;
        }
        case 2: {
            if (cit->second > SCALAR(0)) insert(*cit);
            ++cit;
            break;
        }
        case 3: {
            operator[](it->first) = ((it->second > cit->second) ? (it->second) : (cit->second));
            ++cit;
            ++it;
            break; }
        default:;
        }
    }
    if (cit == source.end())
    {
        for (; it != target.end(); )
            if (!(it->second > SCALAR(0))) erase((it++)->first);
            else ++it;
    }
    if (it == target.end())
    {
        for (; cit != source.end(); ++cit)
            if (cit->second > SCALAR(0)) insert(*cit);
    }
#else
        typename MAP::iterator it(map_begin()), itend(map_end());
        typename MAP::const_iterator cit(rhs.map_begin()), cend(rhs.map_end());
        for (; it != itend && cit != cend;) {
            // c++11 syntax auto
            int c = (it->first < cit->first) ? 1 : (cit->first < it->first) ? 2 : (cit->first == it->first) ? 3 : 4;
            switch (c) {
                case 1: {
                    if (!(it->second > SCALAR(0))) erase(it++);
                    break;
                }
                case 2: {
                    if (cit->second > SCALAR(0)) insert(*cit);
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
        if (cit == cend) {
            for (; it != itend;)
                if (!(it->second > SCALAR(0))) erase(it++);
                else ++it;
        }
        if (it == itend) {
            for (; cit != cend; ++cit)
                if (cit->second > SCALAR(0)) insert(*cit);
        }
#endif // UNORDEREDMAP
        return *this;
    }

    /// Binary version of  operator|=()
    inline __DECLARE_BINARY_OPERATOR(sparse_vector,
                                     |, |=, sparse_vector);

    /// A version of operator+=(rhs.scal_prod(s))
    /// when RHS is a scaled basis vector
    inline sparse_vector &add_scal_prod(const KEY &rhs,
                                        const SCALAR &s)
    {
        // sparse addition
        if (SCALAR(0) == (operator[](rhs) += s)) erase(rhs);
        return *this;
    }

    /// A version of operator+=(rhs.scal_prod(s))
    /// when RHS is a sparse vector scaled on right
    inline sparse_vector &add_scal_prod(const sparse_vector &rhs,
                                        const SCALAR &s)
    {
        if ((s == zero) || rhs.empty())
            return *this;
        if (empty()) {
            *this = rhs;
            return operator*=(s);
        }
        iterator it = begin();
        const_iterator cit;
        for (cit = rhs.begin(); cit != rhs.end(); ++cit) { // Instead of a bare (*this)[cit->first] += cit->second * s;
            it = this->insert(it, std::make_pair(cit->key(),
                                                 zero)); // note this fails if the entry is already there but sets it in any case
            if ((it->value() += cit->value() * s) == zero)
                // erase returns void until c++11
            {
                iterator j(it++);
                erase(j);
            } else ++it;
        }
        return *this;
    }

    /// A version of operator-=(rhs.scal_prod(s))
    /// when RHS is a scaled basis vector
    inline sparse_vector &sub_scal_prod(const KEY &rhs,
                                        const SCALAR &s)
    {
        // sparse addition
        if (SCALAR(0) == (operator[](rhs) -= s)) erase(rhs);
        return *this;
    }

    /// A version of operator-=(rhs.scal_prod(s))
    /// when RHS is a sparse vector scaled on right
    inline sparse_vector &sub_scal_prod(const sparse_vector &rhs,
                                        const SCALAR &s)
    {
        iterator it;
        const_iterator cit;
        if ((s == zero) || rhs.empty())
            return *this;
        if (empty()) {
            *this = rhs;
            return operator*=(-s);
        }
        for (cit = rhs.begin(); cit != rhs.end(); ++cit) { // Instead of a bare (*this)[cit->first] -= cit->second * s;
            it = find(cit->key());
            if (it == end())
                (*this)[cit->key()] = cit->value() * -s;
            else if ((it->value() -= cit->value() * s) == zero)
                erase(it->key());
        }
        return *this;
    }

    /// A fast version of operator+=(rhs.scal_div(s))
    inline sparse_vector &add_scal_div(const sparse_vector &rhs,
                                       const RATIONAL &s)
    {
        iterator it;
        const_iterator cit;
        if (rhs.empty())
            return *this;
        if (empty()) {
            *this = rhs;
            return operator/=(s);
        }
        for (cit = rhs.begin(); cit != rhs.end(); ++cit) { // Instead of a bare (*this)[cit->first] += cit->second / s;
            it = find(cit->key());
            if (it == end())
                (*this)[cit->key()] = cit->value() / s;
            else if ((it->value() += (cit->value() / s)) == zero)
                erase(it->key());
        }
        return *this;
    }

    /// A fast version of operator-=(rhs.scal_div(s))
    inline sparse_vector &sub_scal_div(const sparse_vector &rhs,
                                       const RATIONAL &s)
    {
        iterator it;
        const_iterator cit;
        if (rhs.empty())
            return *this;
        if (empty()) {
            *this = rhs;
            return operator/=(-s);
        }
        for (cit = rhs.begin(); cit != rhs.end(); ++cit) { // Instead of a bare (*this)[cit->first] -= cit->second / s;
            it = find(cit->key());
            if (it == end())
                (*this)[cit->key()] = -cit->value() / s;
            else if ((it->value() -= (cit->value() / s)) == zero)
                erase(it->key());
        }
        return *this;
    }

    inline sparse_vector &add_scal_div(const KEY &rhs, const RATIONAL &s)
    {
        if (zero == (operator[](rhs) += one / s)) erase(rhs);
        return *this;
    }

    inline sparse_vector &sub_scal_div(const KEY &rhs, const RATIONAL &s)
    {
        if (zero == (operator[](rhs) -= one / s)) erase(rhs);
        return *this;
    }

    /// Compares the instance to a sparse_vector.
    bool operator==(const sparse_vector &rhs) const
    {
        if (size() != rhs.size())
            return false;
        const_iterator i, j, jend(rhs.end()), iend(end());
        for (i = begin(); i != iend; ++i) {
            j = rhs.find(i->key());
            if ((j == jend) || (j->value() != i->value()))
                return false;
        }
        return true;
    }

    /// Lexicographically compares the instance to a sparse_vector.
    bool operator<(const sparse_vector &rhs) const
    {
        return std::lexicographical_compare(map_begin(), map_end(), rhs.map_begin(), rhs.map_end());
    }

    /// Boolean negation of operator==()
    bool operator!=(const sparse_vector &rhs) const
    {
        return !operator==(rhs);
    }

    DEG degree() const
    {
        DEG ans(0);
        for (const_iterator it(begin()); it != end(); ++it) {
            ans = std::max(basis.degree(it->key()), ans);
        }
        assert (ans <= BASIS::MAX_DEGREE);
        return ans;
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

    /// Computes the l1 norm of degree d component of a sparse vector with respect to this basis
    inline SCALAR NormL1(const DEG &d) const
    {
        const_iterator i;
        SCALAR ans(zero);
        for (i = begin(); i != end(); ++i) {
            if (d == basis.degree(i->key()))
                ans += abs(i->value());
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

    /// Computes the l-infinity norm of degree d component of a sparse vector with respect to this basis
    inline SCALAR NormLInf(const DEG &d) const
    {
        const_iterator i;
        SCALAR ans(zero);
        for (i = begin(); i != end(); ++i) {
            if (d == basis.degree(i->key()))
                ans = std::max(abs(i->value()), ans);
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

    inline friend std::ostream &operator<<(std::ostream &os,
                                           const sparse_vector &rhs)
    {

        std::pair<BASIS *, KEY> token;
        token.first = &sparse_vector::basis;
        os << '{';
        // create buffer to avoid unnecessary calls to MAP inside loop
#ifndef ORDEREDMAP
        typename std::vector<std::pair < KEY, SCALAR> > ::const_iterator
        cit;
        typename std::vector<std::pair < KEY, SCALAR> >
                buffer(rhs.map_begin(), rhs.map_end());
        std::sort(buffer.begin(), buffer.end(), comp);
        for (cit = buffer.begin(); cit != buffer.end(); ++cit) {
            token.second = cit->first;
            os << ' ' << cit->second << '(' << token << ')';
        }

#else
        const_iterator cit;
        const sparse_vector &buffer = rhs;

        for (cit = rhs.begin(); cit != rhs.end(); ++cit) {
            token.second = cit->key();
            os << ' ' << cit->value() << '(' << token << ')';
        }
#endif // ORDEREDMAP

        os << " }";
        return os;
    }


protected:

    void fill_buffer(
            std::vector<std::pair<KEY, SCALAR> > &buffer
    ) const
    {
        buffer.assign(map_begin(), map_end());
    }

    /// copy the (key, value) elements from rhs to a sorted vector buffer (using the key for sorting)
    /// and construct an increasing vector iterators so that segment [iterators[i-1], iterators[i])
    /// contains keys of degree i; the first begins at [begin(), and the last ends at end), and it can be empty
    void separate_by_degree(
            std::vector<std::pair<KEY, SCALAR> > &buffer,
            const sparse_vector &rhs,
            const size_t DEPTH1,
            std::vector<typename std::vector<std::pair<KEY, SCALAR> >::const_iterator> &iterators
    ) const
    {
        rhs.fill_buffer(buffer);
#ifndef ORDEREDMAP
        std::sort(buffer.begin(), buffer.end(),
                  [](const std::pair<KEY, SCALAR>&lhs, const std::pair<KEY, SCALAR>&rhs)->bool
                  {return lhs.first < rhs.first; }
        );
#endif // ORDEREDMAP

        iterators.assign(DEPTH1 + 1, buffer.end());
        unsigned deg = 0;
        for (typename std::vector<std::pair<KEY, SCALAR> >::const_iterator j0 = buffer.begin();
             j0 != buffer.end();
             j0++) {
            DEG d = basis.degree(j0->first);
            assert(d >= deg && d <= DEPTH1); // order assumed to respect degree
            while (deg < d)
                iterators[deg++] = j0;
            // deg == d
        }
    }

public:
    // Transform methods


    template<typename Vector, typename KeyTransform>
    void triangular_buffered_apply_binary_transform(
            Vector &result,
            const sparse_vector &rhs,
            KeyTransform key_transform,
            const DEG max_depth
    ) const
    {
        // create buffers to avoid unnecessary calls to MAP inside loop
        std::vector<std::pair<KEY, SCALAR> > buffer;
        std::vector<typename std::vector<std::pair<KEY, SCALAR> >::const_iterator>
                iterators;
        separate_by_degree(buffer, rhs, max_depth, iterators);

        typename std::vector<std::pair<KEY, SCALAR> >::const_iterator j, jEnd;
        const_iterator i(begin()), iEnd(end());
        for (; i != iEnd; ++i) {
            const KEY &k = i->key();
            size_t rhdegree = max_depth - basis.degree(k);
            typename std::vector<std::pair<KEY, SCALAR> >::const_iterator &
                    jEnd = iterators[rhdegree];
            for (j = buffer.begin(); j != jEnd; ++j) {
                key_transform(result, i->key(), i->value(), j->first, j->second);
            }
        }
    }

    template<typename Vector, typename KeyTransform, typename IndexTransform>
    void triangular_buffered_apply_binary_transform(
            Vector &result,
            const sparse_vector &rhs,
            KeyTransform key_transform,
            IndexTransform /* index_transform */,
            const DEG max_depth
    ) const
    {
        triangular_buffered_apply_binary_transform(result, rhs, key_transform, max_depth);
    }


    template<typename Vector, typename KeyTransform>
    void square_buffered_apply_binary_transform(
            Vector &result,
            const sparse_vector &rhs,
            KeyTransform key_transform
    ) const
    {
        // create buffer to avoid unnecessary calls to MAP inside loop
        std::vector<std::pair<KEY, SCALAR> > buffer(rhs.map_begin(), rhs.map_end());
        const_iterator i;

        // DEPTH1 == 0
        typename std::vector<std::pair<KEY, SCALAR> >::const_iterator j;
        for (i = begin(); i != end(); ++i) {
            for (j = buffer.begin(); j != buffer.end(); ++j) {
                key_transform(result, i->key(), i->value(), j->first, j->second);
            }
        }


    }

    template<typename Vector, typename KeyTransform, typename IndexTransform>
    void square_buffered_apply_binary_transform(
            Vector &result,
            const sparse_vector &rhs,
            KeyTransform key_transform,
            IndexTransform /* index_transform */
    ) const
    {
        square_buffered_apply_binary_transform(result, rhs, key_transform);
    }


};


} // namespace vectors
} // namespace alg


// Include once wrapper
// DJC_COROPA_LIBALGEBRA_SPARSEVECTORH_SEEN
#endif
//EOF.


