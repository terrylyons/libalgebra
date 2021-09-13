/* *************************************************************

Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurk� and Arend Janssen.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */

//  lie_basis.h

// Include once wrapper
#ifndef DJC_COROPA_LIBALGEBRA_LIEBASISH_SEEN
#define DJC_COROPA_LIBALGEBRA_LIEBASISH_SEEN

#include "libalgebra/basis/basis.h"
#include "utils/integer_maths.h"

namespace dtl {

using alg::integer_maths::mobius;
using alg::integer_maths::power;

/**
 * The size of the Hall set with n letters has number of members at level k
 * given by
 *
 * (1/k) * \sum_{d | k} \mu(d) n^{K/d}
 *
 * where \mu is the Mobius function and the sum is taken over all divisors of d.
 */


template <typename Unsigned>
constexpr std::make_signed_t<Unsigned> hall_set_level_term(Unsigned width, Unsigned level, Unsigned divisor)
{
    using Int = std::make_signed_t<Unsigned>;
    return mobius(divisor)*power(static_cast<Int>(width), level/divisor)/ static_cast<Int>(level);
}

template <typename Unsigned>
constexpr std::make_signed_t<Unsigned> hall_set_level_size_impl(Unsigned width, Unsigned level, Unsigned divisor)
{
    return (
            (level%divisor==0) ? hall_set_level_term(width, level, divisor) : 0
    )+((divisor<level) ? hall_set_level_size_impl(width, level, divisor+1) : 0);
}

template <typename Unsigned>
constexpr Unsigned hall_set_level_size(Unsigned width, Unsigned level)
{
    return static_cast<Unsigned>(hall_set_level_size_impl(width, level, static_cast<Unsigned>(1)));
}

template <typename Unsigned>
constexpr Unsigned hall_set_size_impl(Unsigned width, Unsigned depth, Unsigned level)
{
    return hall_set_level_size(width, level)+
            (level<depth ? hall_set_size_impl(width, depth, level+1) : 0);
}

template <typename Unsigned>
constexpr Unsigned hall_set_size(Unsigned width, Unsigned depth)
{
    return (depth==static_cast<Unsigned>(0))
        ? static_cast<Unsigned>(0)
        : hall_set_size_impl(width, depth, static_cast<Unsigned>(1));
}

template <DIMN Width, DIMN Depth>
struct hall_set_array_helper
{
    static constexpr DIMN value = hall_set_size(Width, Depth);
};

template <DEG NoLetters, DEG MaxDepth> struct hall_set_info
{
    template <DIMN Level>
    using helper = hall_set_array_helper<NoLetters, Level>;

    using holder = typename alg::utils::generate_array<MaxDepth+1, helper>::result;
    static constexpr std::array<DIMN, MaxDepth + 2> degree_sizes = holder::data;
};


} // namespace dtl

/// The Hall Basis class.
/**

   A basis is a finite total ordered set of keys, its cardinal is size() and
   its minimal element is begin(). The successor key of a given key is given
   by nextkey(). The successor of the maximal key is end() and does not belong
   to the basis. The position of a given key in the total order of the basis
   is given by keypos(), and equals 1 for begin(). To each letter corresponds
   a key.

   This class is an ancestor of the lie_basis class, used to implement the lie
   class (Free Lie Algebra) as a particular instance of an algebra class
   (Associative Algebras).

   This class stores a Philip Hall basis associated to a finite number of
   letters. A key is the implementation of a Lie element of this basis. A
   letter is a particular Lie element (or basis element, or key). Each key k
   which does not correspond to a letter has two parents lp and rp and we have
   k = [lp,rp] where [.,.] is the Lie product. A letter, viewed as a key, has
   no parents. More precisely, its parents are invalid keys.

   The basis elements are recursively computed and are enumerated with keys.
   The set of valid keys is essentially an interval of natural integers.

   One can find below a brief Mathematical description of Philip Hall bases
   for the free Lie Algebra. Cf. Reutenauer's book for example, ISBN 0 19
   853679 8.

   Let K be a field with characteristic non equals to 2. In
   newgenesis-libalgebra, this field K corresponds to the type SCA defined in
   libalgebra.h.

   Let M be a finite alphabet {a_1,...,a_n}. We denote by M* the monoid which
   consists in words of letters in M. The product in M* is the concatenation
   and the neutral element is the empty word.

   We consider the free albegra A over (K,M). An element of A is a linear
   combination of elements of M*, with coefficients in K. An element of A is
   an instance of class free_tensor<>, which affects to each element of M* a
   coefficient in K. The element of M* are indexed by tensor_key<>, which
   essentially stores the corresponding word as a std::string.

   We consider also the associated free Lie albegra L, the smallest subalgebra
   of A which contains M and is stable by the Lie product [X,Y] = XY-YX. An
   element of L is an instance of class lie<>. The key used are of type
   lie_key<>, which are actually indexes in a basis of type lie_basis<>.

   The degree of a word w in M is its length. The degree of an element of the
   algebra A is the maximum degree of words with non zero coefficients. The
   degree of [X,Y] is the sum of the degrees of X and Y if X and Y are
   different, and 0 if X = Y.

   Actually, the free Lie algebra L is a graded algebra, with respect to the
   degree (or weight) of Lie products. Philip Hall invented an algorithm for
   computing a basis of the free Lie albegra L. A Hall basis H is a union of
   subsets H_1,...H_i,... of L. By definition, H_1 = M = {a_1,...,a_n} and the
   elements of H_i are of degree i. The set H is totally ordered and more
   over, H_1 < H_2 < ... The Hall basis H can be constructed recursively from
   H_1. This can be done by constructing an array HALLARRAY of elements of the
   form {left, degree , right}. The left and right corresponds to indexes in
   the array for constructing the element by the Lie product, and degree
   corresponds to the degree of this element, which is then the sum of the
   degrees of the two elements pointed by left and right. The order of the
   elements of the array is in one to one correspondance with the order of H.
   The subset H_i is exactly the elements of the form {left, degree , right}
   with degree = i.

   Starting from H1 = {{0, 1, 1},...,{0, 1, n}} which corresponds to the n
   letters, Hi+1 is constructed from H_1, ..., H_i by examining all elements
   of the form {l, i + 1, r} where l < r and l and r are in the union of
   H_1,...,H_i. Such an element is added to the set Hi+1 if and only if the
   right parent of r is <= l.
*/
// Hall basis provides a fully populated (i.e. dense) basis for the lie elements
// has the ability to extend the basis by degrees. This is not protected or
// thread safe. To avoid error in multi-threaded code it is essential that it is
// extended at safe times (e.g. on construction). The code for lie basis has
// been modified to do this and avoid a subtle error.

// It would be worthwhile to write a data driven sparse hall basis

template <DEG N_letters> class hall_basis
{
public:
    /// The default key has value 0, which is an invalid value
    /// and occurs as a parent key of any letter.
    ///
    /// keys can get large - but in the dense case this is not likely
    /// Make a choice for the length of a key in 64 bit.
    // typedef DEG KEY; // unsigned int
    typedef LET KEY; // size_t
    /// The parents of a key are a pair of prior keys. Invalid 0 keys for letters.
    typedef std::pair<KEY, KEY> PARENT;
    /// The number of letters in alphabet
    enum { n_letters = N_letters };

protected:
    /// Parents, indexed by keys.
    std::vector<PARENT> hall_set;

public:
    /// Reverse map from parents to keys.
    struct REVERSE_MAP : private std::map<PARENT, KEY>
    {
        using std::map<PARENT, KEY>::operator[];
        using std::map<PARENT, KEY>::find;
        using std::map<PARENT, KEY>::end;
    };

protected:
    REVERSE_MAP reverse_map;
    /// Degrees, indexed by keys.
    // static struct DEGREE : private
    struct DEGREE : private std::vector<DEG>
    {
        using std::vector<DEG>::operator[];
        using std::vector<DEG>::push_back;
        using std::vector<DEG>::begin;
        using std::vector<DEG>::resize;
    } degrees;
    /// Index to basis pointing at the first and one past the last of a given
    /// degree
    std::vector<std::pair<size_t, size_t>> hall_set_degree_ranges;
    /// Letters, indexed by their keys.
    // std::string letters;
    std::vector<LET> letters;
    /// Maps letters to keys.
    std::map<LET, KEY> ltk;
    /// Current degree, always > 0 for any valid class instance.
    DEG curr_degree;

public:
    /// Constructs the basis with a given number of letters.
    hall_basis() : curr_degree(0)
    {
        // We put something at position 0 to make indexing more logical

        degrees.push_back(0);
        hall_set.push_back(PARENT(0, 0));
        hall_set_degree_ranges.push_back(std::pair<DIMN, DIMN>(0, hall_set.size()));

        for (LET c = 1; c <= n_letters; ++c) {
            letters.push_back(c);
        } //+= (char) c;

        // We add the letters to the basis, starting from position 1.

        for (LET i = 1; i <= letters.size(); ++i) {
            PARENT parents(0, i);
            hall_set.push_back(parents); // at [i]
            reverse_map[parents] = hall_set.size() - 1;
            degrees.push_back(1); // at [i]
            ltk[letters[i - 1]] = (LET) i;
        }
        std::pair<size_t, size_t> range;
        range.first = hall_set_degree_ranges[curr_degree].second;
        range.second = hall_set.size();
        hall_set_degree_ranges.push_back(range);
        ++curr_degree;
        // To construct the rest of the basis now, call growup(max_degree) here.
    }
    /// Constructs the basis up to a desired degree.
    /**
    For performance reasons, max_degree is not checked. So be careful.
    */
    inline void growup(DEG desired_degree)
    {
        for (DEG d = curr_degree + 1; d <= desired_degree; ++d) {
            for (DEG e = 1; 2 * e <= d; ++e) {
                LET i_lower = hall_set_degree_ranges[e].first;
                LET i_upper = hall_set_degree_ranges[e].second;
                LET j_lower = hall_set_degree_ranges[d - e].first;
                LET j_upper = hall_set_degree_ranges[d - e].second;
                for (LET i = i_lower; i < i_upper; ++i) {
                    for (LET j = std::max(j_lower, i + 1); j < j_upper; ++j) {
                        if (hall_set[j].first <= i) {
                            PARENT parents(i, j);
                            hall_set.push_back(parents);
                            degrees.push_back(d);
                            assert(((degrees[i] + degrees[j]) == degrees[hall_set.size() - 1]));
                            reverse_map[parents] = hall_set.size() - 1;
                        }
                    }
                }
            }
            std::pair<size_t, size_t> range;
            range.first = hall_set_degree_ranges[curr_degree].second;
            range.second = hall_set.size();
            hall_set_degree_ranges.push_back(range);
            ++curr_degree;
        }
    }

    /// Returns the degree (ie. weight) of a Lie key.
    inline DEG degree(const KEY &k) const
    {
        assert(k <= size());
        return degrees[k];
    }

    /// Returns the key corresponding to a letter.
    inline KEY keyofletter(LET letter) const { return ltk.find(letter)->second; }

    /// Returns the left parent of a key.
    inline KEY lparent(const KEY &k) const { return hall_set[k].first; }

    /// Returns the right parent of a key.
    inline KEY rparent(const KEY &k) const { return hall_set[k].second; }

    /// Tells if a key corresponds to a letter.
    inline bool letter(const KEY &k) const
    {
        return ((k > 0) && (k <= letters.size()));
    }

    /// Returns the letter of a key corresponding to a letter.
    inline LET getletter(const KEY &k) const { return letters[k - 1]; }

    /// Returns the value of the smallest key in the basis.
    inline KEY begin(void) const { return 1; }

    /// Returns the key next the biggest key of the basis.
    inline KEY end(void) const { return 0; }

    /// Returns the key next a given key in the basis. No implicit growup made.
    inline KEY nextkey(const KEY &k) const
    {
        if (k < (hall_set.size() - 1)) {
            return (k + 1);
        } else {
            return 0;
        }
    }

    /// Returns the position of a key in the basis total order.
    inline DEG keypos(const KEY &k) const { return k; }

    /// Returns the size of the basis.
    inline KEY size(void) const { return (hall_set.size() - 1); }

    /// Outputs the Hall basis to an std::ostream.
    inline friend std::ostream &operator<<(std::ostream &os, const hall_basis &b)
    {
        for (KEY k = b.begin(); k != b.end(); k = b.nextkey(k)) {
            os << b.key2string(k) << ' ';
        }
        return os;
    }

    // inline const std::string& key2string(const KEY& k) const
    // BUG//TJL//24/08/2012 - returned reference invalidated if vector grows!!
    // BUG//TJL//25/08/2012 - not templated but has static member so this is
    // shared across all dimensions regardless of no letters etc!!
private:
    mutable std::vector<std::string> table; // move to instance per class
public:
    // ShortFix return a value not a reference
    // TODO check performance of fix 24/08/2012
    /// Converts a key to an std::string of letters.

    inline const std::string key2string(const KEY &k) const
    {
        static boost::recursive_mutex table_access;
        //// get exclusive recursive access for the thread
        boost::lock_guard<boost::recursive_mutex> lock(table_access);

        // BUG//TJL//09/04/2017 - non static member added to class but not commented
        // out here!!
        //		static std::vector<std::string> table;

        if (k > table.size()) {
            for (KEY i = (KEY) table.size() + 1; i <= k; ++i) {
                table.push_back(_key2string(i));
            }
        }
        return table[k - 1];
    }

private:
    /// Recursively constructs the string associated to the Lie key k.
    std::string _key2string(const KEY &k) const
    {
        std::ostringstream oss;
        if (k > 0) {
            if (letter(k)) {
                oss << getletter(k);
            } else {
                oss << '[';
                oss << key2string(lparent(k));
                oss << ',';
                oss << key2string(rparent(k));
                oss << ']';
            }
        }
        return oss.str();
    }


public:

    struct no_caching_tag {};

    template <DEG CacheDepth>
    struct lazy_cache_tag {};

    template <DEG CacheDepth>
    struct lookup_table_tag {};

    template <typename Function, typename BinOp, typename Tag>
    class extended_function
    {
    public:
        typedef Function function_type;
        typedef BinOp binary_operation_type;
        typedef Tag tag_type;
        typedef decltype(std::declval<function_type>()(std::declval<KEY>())) output_type;

        using table_t = std::unordered_map<KEY, output_type>;

        extended_function(hall_basis& hs) : m_hall_set(hs), m_tag(), m_fn(), m_op()
        {}

        extended_function(Function fn, BinOp op, hall_basis& hs) : m_hall_set(hs), m_tag(), m_fn(fn), m_op(op)
        {}

    private:

        output_type eval_impl(const KEY& k) const
        {
            if (m_hall_set.letter(k)) {
                return m_fn(m_hall_set.getletter(k));
            } else {
                return m_op(
                        operator()(m_hall_set.lparent(k)),
                        operator()(m_hall_set.rparent(k))
                        );
            }
        }

        output_type eval(const KEY& k, no_caching_tag) const
        {
            return eval_impl(k);
        }

        template <DEG CacheDepth>
        output_type eval(const KEY& k, lazy_cache_tag<CacheDepth>) const
        {
            static boost::recursive_mutex table_lock;
            static table_t table;

            if (m_hall_set.degree(k) <= CacheDepth) {
                boost::lock_guard<boost::recursive_mutex> access(table_lock);

                typename table_t::iterator it = table.find(k);
                if (it != table.end()) {
                    return it->second;
                }

                return table[k] = eval_impl(k);
            } else {
                return eval_impl(k);
            }
        }

        output_type eval(const KEY& k, lazy_cache_tag<0>) const
        {
            static boost::recursive_mutex table_lock;
            static table_t table;

            boost::lock_guard<boost::recursive_mutex> access(table_lock);

            typename table_t::iterator it = table.find(k);
            if (it != table.end()) {
                return it->second;
            }

            return table[k] = eval_impl(k);
        }

        template <DEG Depth>
        table_t fill_table() const
        {
            table_t result;

            KEY k = 1;
            for (; m_hall_set.degree(k)==1; ++k) {
                result[k] = m_fn(m_hall_set.getletter(k));
            }

            for (; m_hall_set.degree(k) <= Depth; ++k) {
                result[k] = m_op(result[m_hall_set.lparent(k)], result[m_hall_set.rparent(k)]);
            }

            return result;
        }

        template <DEG CacheDepth>
        output_type eval(const KEY& k, lookup_table_tag<CacheDepth>) const
        {
            static table_t table = fill_table<CacheDepth>();

            if (m_hall_set.degree(k) <= CacheDepth) {
                return table[k];
            } else {
                return eval_impl(k);
            }
        }


    public:

        output_type operator()(const KEY& k) const
        {
            return eval(k, m_tag);
        }

    private:
        hall_basis& m_hall_set;
        function_type m_fn;
        binary_operation_type m_op;
        tag_type m_tag;
    };

};
//// if degree is static
// template<DEG n_letters>
// typename hall_basis<n_letters>::DEGREE hall_basis<n_letters>::degrees;

/// The Lie basis class.
/**
 This is the basis used to implement the lie class as a specialisation of
 the algebra class. In the current implementation, the Lie basis class is a
 wrapper for the hall_basis class, with a prod() member function.
*/

template <DEG n_letters, DEG max_degree>
class lie_basis : protected hall_basis<n_letters>, dtl::hall_set_info<n_letters, max_degree>
{
    typedef dtl::hall_set_info<n_letters, max_degree> SIZE_INFO;

public:
    /// Import of the KEY type.
    typedef typename hall_basis<n_letters>::KEY KEY;
    /// Import of the PARENT type.
    typedef typename hall_basis<n_letters>::PARENT PARENT;
    /// Import member functions.
    using hall_basis<n_letters>::letter;
    using hall_basis<n_letters>::getletter;
    using hall_basis<n_letters>::lparent;
    using hall_basis<n_letters>::rparent;
    // using hall_basis<n_letters>::degrees;
    using hall_basis<n_letters>::growup;
    using hall_basis<n_letters>::reverse_map;
    using hall_basis<n_letters>::degree;
    using hall_basis<n_letters>::keyofletter;
    using hall_basis<n_letters>::begin;
    using hall_basis<n_letters>::end;
    using hall_basis<n_letters>::nextkey;
    using hall_basis<n_letters>::key2string;
    using hall_basis<n_letters>::size;
    using hall_basis<n_letters>::hall_set;

    using typename hall_basis<n_letters>::REVERSE_MAP;

public:
    typedef basis::with_degree<max_degree> degree_tag;
    typedef basis::ordered<std::less<KEY>> ordering_tag;

public:
    /// Constructs the basis for a finite number of letters.
    lie_basis(void) : hall_basis<n_letters>()
    {
        // bug: tjl : 08 04 2017 without the following line the basis would not
        // remain const and sharing it between threads would cause errors
        hall_basis<n_letters>::growup(max_degree);
    }

    /// Returns the product of two key.
    /**
    Returns the LIE instance corresponding to the product of the two basis
    elements k1 and k2. For performance reasons, the basis is enriched/grown
    dynamically and the already computed products are stored in a static
    multiplication table to speed up further calculations. This function
    returns a constant reference to the suitable table element. See Ch4 p93 and 94
    Reutenauer
    */

    /// Replaces letters by lie<> instances in a lie<> instance.
    /**
    Replaces the occurrences of s letters in the expression of k by the lie<>
    elements in v, and returns the recursively expanded result. The already
    computed replacements are stored in table.
    */
    template <typename Lie>
    Lie replace(const KEY &k, const std::vector<LET> &s, const std::vector<const Lie *> &v, std::map<KEY, Lie> &table)
    {
        typename std::map<KEY, Lie>::iterator it;
        it = table.find(k);
        if (it != table.end()) {
            return it->second;
        } else {
            if (letter(k)) {
                typename std::vector<LET>::size_type i;
                for (i = 0; i < s.size(); ++i) {
                    if (s[i] == getletter(k)) {
                        return table[k] = *(v[i]);
                    }
                }
                return (table[k] = (Lie) k);
            } else {
                return (table[k] = replace(lparent(k), s, v, table) * replace(rparent(k), s, v, table));
            }
        }
    }

    /// Outputs the lie basis to an std::ostream.
    inline friend std::ostream &operator<<(std::ostream &os, const lie_basis &b)
    {
        return os << (const hall_basis<n_letters> &) b;
    }

    /// Outupts an std::pair<lie_basis*, KEY> to an std::ostream.
    inline friend std::ostream &operator<<(std::ostream &os, const std::pair<lie_basis *, KEY> &t)
    {
        return os << t.first->key2string(t.second);
    }

public:
    // index_to_key and friends

    /// Convert a key to index in basis order
    static DIMN key_to_index(const KEY k) { return DIMN(k - 1); }

    /// Convert an index to key
    static KEY index_to_key(const DIMN idx) { return KEY(idx + 1); }

    static DIMN start_of_degree(const DEG d)
    {
        assert(d <= max_degree + 1);
        return (d == 0) ? 0 : SIZE_INFO::degree_sizes[d-1];
    }
};

// Include once wrapper
#endif // DJC_COROPA_LIBALGEBRA_LIEBASISH_SEEN

// EOF.
