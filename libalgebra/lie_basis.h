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
 * (1/k) * \\sum_{d | k} \\mu(d) n^{K/d}
 *
 * where \\mu is the Mobius function and the sum is taken over all divisors of d.
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

    using hall_set_type = hall_set<N_letters>;

    using hall_set_type::extended_function;

    using typename hall_set_type::parent_type;

private:

    static DEG letter_degree(LET) noexcept
    {
        return 1;
    }

    struct depth_predicate
    {
        bool operator()(const KEY& key) const noexcept
        {
            auto bound = hall_set_type::start_of_degree(MaxDepth+1);
            return key <= bound;
        }
    };

    using degree_func_type = typename hall_set_type::template extended_function<
            decltype(&letter_degree),
            std::plus<DEG>,
            lookup_table_tag<depth_predicate>
            >;

public:

    const degree_func_type degree;

    /// Constructs the basis with a given number of letters.
    hall_basis() : hall_set_type(MaxDepth), degree(*this, &letter_degree, std::plus<DEG>())
    {
        // Instantiate the degree function table.
        degree(KEY(size() - 1));
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

    /// Returns the key next a given key in the basis. No implicit growup made.
    inline KEY nextkey(const KEY &k) const
    {
        if (k < (hall_set_type::start_of_degree(MaxDepth+1) - 1)) {
            return (k + 1);
        } else {
            return 0;
        }
    }

    friend std::ostream& operator<<(std::ostream& os, const hall_basis& b)
    {
        for (KEY k = b.begin(); k!=b.end(); k = b.next_key(k)) {
            os << b.key2string(k) << ' ';
        }
        return os;
    }

    DIMN size() const
    {
        return start_of_degree(MaxDepth+1);
    }


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

    static DIMN start_of_degree(const DEG d)
    {
        assert(d <= max_degree + 1);
        return (d == 0) ? 0 : SIZE_INFO::degree_sizes[d-1];
    }
};

// Include once wrapper
#endif // DJC_COROPA_LIBALGEBRA_LIEBASISH_SEEN

// EOF.
