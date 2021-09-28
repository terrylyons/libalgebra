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

    using hall_set_type::parent_type;

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
    }

    using hall_set_type::keyofletter;
    using hall_set_type::lparent;
    using hall_set_type::rparent;
    using hall_set_type::letter;
    using hall_set_type::getletter;
    using hall_set_type::key2string;
    using hall_set_type::size;
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
        return os << static_cast<const hall_set_type&>(b);
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
