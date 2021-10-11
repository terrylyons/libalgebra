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
#include "libalgebra/utils/integer_maths.h"

namespace dtl {

using alg::integer_maths::divisor_calc;
using alg::integer_maths::mobius_func;

typedef long long Long;

/**
 * The size of the Hall set with n letters has number of members at level k
 * given by
 *
 * (1/k) * \\sum_{d | k} \\mu(d) n^{K/d}
 *
 * where \\mu is the Mobius function and the sum is taken over all divisors of d.
 */

/*
Struct for computing the size of a layer in the Hall set.

We recurse down through divisors, adding on the corresponding
term from the formula for the size of the layer. This quantity
is not normalised.
*/
template <DEG NoLetters, DEG Level, DEG Divisor = 1, DEG Remainder = (Level % Divisor)> struct hall_set_level_size
{
    //   static constexpr const Long value = hall_set_level_size<NoLetters, Level,
    //   Divisor + 1>::value;
    enum { value = hall_set_level_size<NoLetters, Level, Divisor + 1>::value };
};

/*
Specialisation for divisors: multiply mobius by exponent and add to next
divisor term value.
*/
template <DEG NoLetters, DEG Level, DEG Divisor> struct hall_set_level_size<NoLetters, Level, Divisor, 0>
{
    /*
    static constexpr const Long value = mobius_func<Divisor>::value
                              * ConstPower<NoLetters, Level / Divisor>::ans
                              + hall_set_level_size<NoLetters, Level, Divisor +
    1>::value;
                              */
    enum
    {
        value = mobius_func<Divisor>::value * ConstPower<NoLetters, Level / Divisor>::ans +
                hall_set_level_size<NoLetters, Level, Divisor + 1>::value
    };
};

/*
Specialisation for divisor = number (terminal case).
*/
template <DEG NoLetters, DEG Level> struct hall_set_level_size<NoLetters, Level, Level, 0>
{
    // static constexpr const Long value = mobius_func<Level>::value * NoLetters;
    enum { value = mobius_func<Level>::value * NoLetters };
};

/*
Specialisation for 1
*/
template <DEG NoLetters> struct hall_set_level_size<NoLetters, 1, 1, 0>
{
    // static constexpr const Long value = NoLetters;
    enum { value = NoLetters };
};

/*
Specialisation for 0.
*/
template <DEG NoLetters, DEG Divisor> struct hall_set_level_size<NoLetters, 0, Divisor, 0>
{
    // static constexpr const Long value = NoLetters;
    enum { value = NoLetters };
};

template <DEG NoLetters, DEG MaxLevel> struct hall_set_size
{
    // static const Long value;

    // static constexpr const Long value = ((hall_set_level_size<NoLetters,
    // MaxLevel - 1>::value / (MaxLevel - 1))
    //     //    + hall_set_size<NoLetters, MaxLevel - 1>::value);
    enum
    {
        value = ((hall_set_level_size<NoLetters, MaxLevel - 1>::value / (MaxLevel - 1)) +
                 hall_set_size<NoLetters, MaxLevel - 1>::value)
    };
};

template <DEG NoLetters> struct hall_set_size<NoLetters, 1>
{
    // static constexpr const Long value = 0;
    enum { value = 0 };
};

template <DEG NoLetters> struct hall_set_size<NoLetters, 0>
{
    // static constexpr const Long value = 0;
    enum { value = 0 };
};

/*
template < DEG NoLetters, DEG MaxLevel >
const Long hall_set_size < NoLetters, MaxLevel >::value =
        ((hall_set_level_size<NoLetters, MaxLevel - 1>::value / (MaxLevel - 1))
          + hall_set_size<NoLetters, MaxLevel - 1>::value);

template < DEG NoLetters >
const Long hall_set_size < NoLetters, 1 >::value = 0;

template < DEG NoLetters >
const Long hall_set_size < NoLetters, 0 >::value = 0;
*/

template <DEG NoLetters, DEG MaxDepth> struct hall_set_info
{
    static const LIBALGEBRA_STATIC_ARRAY_TYPE<DIMN, MaxDepth + 2> degree_sizes;
};

template <DEG NoLetters, DEG MaxDepth> LIBALGEBRA_STATIC_ARRAY_TYPE<DIMN, MaxDepth + 2> populate_hall_set_size_array()
{
    LIBALGEBRA_STATIC_ARRAY_TYPE<DIMN, MaxDepth + 2> tmp;
    utils::populate_array<hall_set_size, NoLetters, MaxDepth + 1>::fill(tmp);
    return tmp;
}

template <DEG NoLetters, DEG MaxDepth> const LIBALGEBRA_STATIC_ARRAY_TYPE<DIMN, MaxDepth + 2>
        hall_set_info<NoLetters, MaxDepth>::degree_sizes = populate_hall_set_size_array<NoLetters, MaxDepth>();

} // namespace dtl

// Hall basis provides a fully populated (i.e. dense) basis for the lie elements
// has the ability to extend the basis by degrees. This is not protected or
// thread safe. To avoid error in multi-threaded code it is essential that it is
// extended at safe times (e.g. on construction). The code for lie basis has
// been modified to do this and avoid a subtle error.

// It would be worthwhile to write a data driven sparse hall basis
/**
 * @brief The Hall basis class
 *
 * The Hall basis is a specific realisation of a Hall set, with N_letters
 * and including elements up to and including MaxDepth. This class is a
 * thin wrapper around the more general hall_set class and exposes some
 * functionality of this underlying class. However, one cannot access
 * information about the underlying hall_set beyond the MaxDepth set in
 * the template parameter, even if the underlying hall_set object does
 * include these data.
 *
 * @tparam N_letters Size of alphabet
 * @tparam MaxDepth Maximum degree of elements to consider
 */
template <DEG N_letters, DEG MaxDepth>
class hall_basis : private hall_set<N_letters>
{
    /// Base class - privately derived, so keep private
    using hall_set_type = hall_set<N_letters>;
public:

    using typename hall_set_type::letter_type;

    /// The default key has value 0, which is an invalid value
    /// and occurs as a parent key of any letter.
    ///
    /// keys can get large - but in the dense case this is not likely
    /// Make a choice for the length of a key in 64 bit.
    using typename hall_set_type::key_type;

    /// The parents of a key are a pair of prior keys. Invalid 0 keys for letters.
    using typename hall_set_type::parent_type;

    /// The number of letters in alphabet
    enum { n_letters = N_letters };

    /// Bring up the extended_function class - see also extend_function method
    template <typename Function, typename BinOp, typename Tag>
    using extended_function = typename hall_set_type::template
            extended_function<Function, BinOp, Tag>;

    /// For external compatibility, define KEY = key_type
    typedef LET KEY; // size_t
    /// For external compatibility, define PARENT = parent_type
    typedef parent_type PARENT;

private:

    /// Key level function for getting the degree
    static DEG letter_degree(letter_type) noexcept
    {
        return 1;
    }

    /// Predicate for forming the lookup table for degree
    template <DEG D=MaxDepth>
    struct depth_predicate
    {
        bool operator()(const key_type& key) const noexcept
        {
            auto bound = hall_set_type::start_of_degree(D+1);
            return key <= bound;
        }
    };

    /// Extended function for computing degree
    using degree_func_type = typename hall_set_type::template extended_function<
            decltype(&letter_degree),
            std::plus<DEG>,
            lookup_table_tag<depth_predicate<MaxDepth> >
            >;

public:

    /// Function returning the degree of a given key
    const degree_func_type degree;

    /// Constructs the basis with a given number of letters.
    hall_basis() : hall_set_type(MaxDepth), degree(*this, &letter_degree, std::plus<DEG>())
    {
        // Instantiate the degree function table.
        degree(letter_type(1));
    }

    using hall_set_type::keyofletter;
    using hall_set_type::lparent;
    using hall_set_type::rparent;
    using hall_set_type::letter;
    using hall_set_type::getletter;
    using hall_set_type::key2string;
    using hall_set_type::operator[];
    using hall_set_type::letters;
    using hall_set_type::hall_set_degree_ranges;
    using hall_set_type::l2k;
    using hall_set_type::find;
    using hall_set_type::reverse_map_end;
    using hall_set_type::parents_begin;
    using hall_set_type::parents_end;
    using hall_set_type::start_of_degree;


    /// Returns the value of the smallest key in the basis.
    using hall_set_type::begin;

    /// Returns the key next the biggest key of the basis.
    using hall_set_type::end;


    /**
     * @brief Returns the key next a given key in the basis.
     * @param k current key
     * @return key that immediately follows k in the Hall set order
     */
    inline KEY nextkey(const KEY &k) const
    {
        if (k < (hall_set_type::start_of_degree(MaxDepth+1) - 1)) {
            return (k + 1);
        } else {
            return 0;
        }
    }

    /// Print the basis to an output stream
    friend std::ostream& operator<<(std::ostream& os, const hall_basis& b)
    {
        for (KEY k = b.begin(); k!=b.end(); k = b.next_key(k)) {
            os << b.key2string(k) << ' ';
        }
        return os;
    }

    /// The size of the Hall set up to degree MaxDepth
    DIMN size() const
    {
        return start_of_degree(MaxDepth+1);
    }

    /**
     * @brief Extend a function on letters to a function on the Hall set
     *
     * Given a function f from letters into some set M with a binary operation
     * * we can define an extension of f to the Hall set H by taking
     * f([k1, k2]) = f(k1) * f(k2) for keys k1 and k2. By the definition of the
     * Hall set, this recursion will terminate when k1 and k2 are letters.\n\n
     *
     * An easy example of this process is the degree function d, which maps
     * letters to the integer 1. The binary operation here is +, and the
     * extension d([k1, k2]) = d(k1) + d(k2) defines the extension to the
     * whole Hall set.\n\n
     *
     * This function is a convenience for constructing such an extension. Given
     * a function-like object Function, a binary operation BinOp, and a tag to
     * indicate the type of caching to use, we construct an extended function
     * object that behaves like the Function on letters.\n\n
     *
     * The tag object can be one three different options: no_caching_tag for no
     * caching; lazy_caching_tag for lazy (cache on first compute) caching; and
     * lookup_table_tag for precomputed lookup table caching.
     *
     * @tparam Function Function-like type on letters
     * @tparam BinOp Binary operation type
     * @tparam Tag Caching indicator tag
     * @param fn Function-like instance
     * @param op Binary operation object
     * @return extended_function instance
     */
    template <typename Tag, typename Function, typename BinOp>
    typename hall_set_type::template extended_function<Function, BinOp, Tag>
    extend_function(Function fn, BinOp op) const
    {
        using ext_fn_t = typename hall_set_type::template extended_function<Function, BinOp, Tag>;
        return ext_fn_t(*this, fn, op);
    }

    /**
     * @brief Extend a function on letters to a function on the Hall set with lazy caching based on depth
     *
     * @tparam Depth Maximum degree of keys to include in cache
     * @tparam Function Function-like type on letters
     * @tparam BinOp Binary operation type
     * @param fn Function-like instance
     * @param op Binary operation object
     * @return extended_function instance
     */
    template <DEG Depth, typename Function, typename BinOp>
    typename hall_set_type::template
            extended_function<Function, BinOp, lazy_cache_tag<depth_predicate<Depth> > >
    extend_function_lazy_cache(Function fn, BinOp op) const
    {
        using tag_type = lazy_cache_tag<depth_predicate<Depth> >;
        using ext_fn_t = typename hall_set_type::template extended_function<Function, BinOp, tag_type>;
        return ext_fn_t(*this, fn, op);
    }

    template <typename ExtendedFunction>
    ExtendedFunction extend_function() const
    {
        return ExtendedFunction(*this);
    }


};

/**
* @brief The Lie basis class.
*
* This is the basis used to implement the lie class as a specialisation of
* the algebra class. In the current implementation, the Lie basis class is a
* wrapper for the hall_basis class, with a prod() member function.
*
* @tparam n_letters Size of alphabet
* @tparam max_degree maximum degree of elements
*/
template <DEG n_letters, DEG max_degree>
class lie_basis : protected hall_basis<n_letters, max_degree>, dtl::hall_set_info<n_letters, max_degree>
{
    typedef dtl::hall_set_info<n_letters, max_degree> SIZE_INFO;

    using hall_basis_type = hall_basis<n_letters, max_degree>;
public:
    /// Import of the KEY type.
    typedef typename hall_basis_type::KEY KEY;
    /// Import of the PARENT type.
    typedef typename hall_basis_type::PARENT PARENT;
    /// Import member functions.
    using hall_basis_type::letter;
    using hall_basis_type::getletter;
    using hall_basis_type::lparent;
    using hall_basis_type::rparent;


    using hall_basis_type::degree;
    using hall_basis_type::keyofletter;
    using hall_basis_type::begin;
    using hall_basis_type::end;
    using hall_basis_type::nextkey;
    using hall_basis_type::key2string;
    using hall_basis_type::size;

    using hall_basis_type::operator[];
    using hall_basis_type::find;
    using hall_basis_type::reverse_map_end;
    using hall_basis_type::parents_begin;
    using hall_basis_type::parents_end;


    /// Bring up the extended_function class - see also extend_function method
    template<typename Function, typename BinOp, typename Tag>
    using extended_function = typename hall_basis_type::template
    extended_function<Function, BinOp, Tag>;

    using hall_basis_type::extend_function;
    using hall_basis_type::extend_function_lazy_cache;

public:
    typedef basis::with_degree<max_degree> degree_tag;
    typedef basis::ordered<std::less<KEY>> ordering_tag;

public:
    /// Constructs the basis for a finite number of letters.
    lie_basis(void) : hall_basis_type()
    {
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
        return os << (const hall_basis_type &) b;
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

    /// Get the index at which elements of given degree start
    static DIMN start_of_degree(const DEG d)
    {
        assert(d <= max_degree + 1);
        auto& r = SIZE_INFO::degree_sizes;
        auto p = hall_basis_type::start_of_degree(d);
        assert(p == r[d]);
        return SIZE_INFO::degree_sizes[d];
    }
};

// Include once wrapper
#endif // DJC_COROPA_LIBALGEBRA_LIEBASISH_SEEN

// EOF.
