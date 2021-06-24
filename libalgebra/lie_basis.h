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
 * The size of the Hall set with n letters has number of members at level k given by
 *
 * (1/k) * \sum_{d | k} \mu(d) n^{K/d}
 *
 * where \mu is the Mobius function and the sum is taken over all divisors of d.
 */

/*
Struct for computing the size of a layer in the Hall set.

We recurse down through divisors, adding on the corresponding
term from the formula for the size of the layer. This quantity
is not normalised.
*/
template<DEG NoLetters, DEG Level, DEG Divisor = 1, DEG Remainder = (Level % Divisor)>
struct hall_set_level_size
{
 //   static constexpr const Long value = hall_set_level_size<NoLetters, Level, Divisor + 1>::value;
    enum {
        value = hall_set_level_size<NoLetters, Level, Divisor + 1>::value
    };
};

/*
Specialisation for divisors: multiply mobius by exponent and add to next
divisor term value.
*/
template<DEG NoLetters, DEG Level, DEG Divisor>
struct hall_set_level_size<NoLetters, Level, Divisor, 0>
{
    /*
    static constexpr const Long value = mobius_func<Divisor>::value
                              * ConstPower<NoLetters, Level / Divisor>::ans
                              + hall_set_level_size<NoLetters, Level, Divisor + 1>::value;
                              */
    enum {
        value = mobius_func<Divisor>::value
                * ConstPower<NoLetters, Level / Divisor>::ans
                + hall_set_level_size<NoLetters, Level, Divisor + 1>::value
    };
};

/*
Specialisation for divisor = number (terminal case).
*/
template<DEG NoLetters, DEG Level>
struct hall_set_level_size<NoLetters, Level, Level, 0>
{
    //static constexpr const Long value = mobius_func<Level>::value * NoLetters;
    enum {
        value = mobius_func<Level>::value * NoLetters
    };
};

/*
Specialisation for 1
*/
template<DEG NoLetters>
struct hall_set_level_size<NoLetters, 1, 1, 0>
{
    //static constexpr const Long value = NoLetters;
    enum {
        value = NoLetters
    };
};



/*
Specialisation for 0.
*/
template<DEG NoLetters, DEG Divisor>
struct hall_set_level_size<NoLetters, 0, Divisor, 0>
{
   // static constexpr const Long value = NoLetters;
   enum {
       value = NoLetters
   };
};



template<DEG NoLetters, DEG MaxLevel>
struct hall_set_size
{
    //static const Long value;

    //static constexpr const Long value = ((hall_set_level_size<NoLetters, MaxLevel - 1>::value / (MaxLevel - 1))
    //    //    + hall_set_size<NoLetters, MaxLevel - 1>::value);
    enum {
        value = ((hall_set_level_size<NoLetters, MaxLevel - 1>::value / (MaxLevel - 1))
            + hall_set_size<NoLetters, MaxLevel - 1>::value)
    };

};

template<DEG NoLetters>
struct hall_set_size<NoLetters, 1>
{
    //static constexpr const Long value = 0;
    enum {
        value = 0
    };
};

template<DEG NoLetters>
struct hall_set_size<NoLetters, 0>
{
    //static constexpr const Long value = 0;
    enum {
        value = 0
    };
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

template<DEG NoLetters, DEG MaxDepth>
struct hall_set_info
{
    static const LIBALGEBRA_STATIC_ARRAY_TYPE<DIMN, MaxDepth + 2> degree_sizes;
};


template<DEG NoLetters, DEG MaxDepth>
LIBALGEBRA_STATIC_ARRAY_TYPE<DIMN, MaxDepth + 2>
populate_hall_set_size_array()
{
    LIBALGEBRA_STATIC_ARRAY_TYPE<DIMN, MaxDepth + 2> tmp;
    utils::populate_array<hall_set_size, NoLetters, MaxDepth+1>::fill(tmp);
    return tmp;
}

template<DEG NoLetters, DEG MaxDepth>
const LIBALGEBRA_STATIC_ARRAY_TYPE<DIMN, MaxDepth + 2> hall_set_info<NoLetters, MaxDepth>::degree_sizes
        = populate_hall_set_size_array<NoLetters, MaxDepth>();

}





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
// has the ability to extend the basis by degrees. This is not protected or thread safe.
// To avoid error in multi-threaded code it is essential that it is extended at safe times (e.g.
// on construction). The code for lie basis has been modified to do this and avoid a subtle error.

// It would be worthwhile to write a data driven sparse hall basis

template<DEG N_letters>
class hall_basis
{
public:
    /// The default key has value 0, which is an invalid value
    /// and occurs as a parent key of any letter.
    ///
    /// keys can get large - but in the dense case this is not likely
    /// Make a choice for the length of a key in 64 bit.
    //typedef DEG KEY; // unsigned int
    typedef LET KEY; // size_t
    /// The parents of a key are a pair of prior keys. Invalid 0 keys for letters.
    typedef std::pair<KEY, KEY> PARENT;
    /// The number of letters in alphabet
    enum
    {
        n_letters = N_letters
    };
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
    //static struct DEGREE : private
    struct DEGREE : private std::vector<DEG>
    {
        using std::vector<DEG>::operator[];
        using std::vector<DEG>::push_back;
        using std::vector<DEG>::begin;
        using std::vector<DEG>::resize;
    } degrees;
    /// Index to basis pointing at the first and one past the last of a given degree
    std::vector<std::pair<size_t, size_t> > hall_set_degree_ranges;
    /// Letters, indexed by their keys.
    //std::string letters;
    std::vector<LET> letters;
    /// Maps letters to keys.
    std::map<LET, KEY> ltk;
    /// Current degree, always > 0 for any valid class instance.
    DEG curr_degree;
public:
    /// Constructs the basis with a given number of letters.
    hall_basis()
            : curr_degree(0)
    {
        // We put something at position 0 to make indexing more logical

        degrees.push_back(0);
        hall_set.push_back(PARENT(0, 0));
        hall_set_degree_ranges.push_back(std::pair<DIMN, DIMN>(0, hall_set.size()));

        for (LET c = 1; c <= n_letters; ++c)
            letters.push_back(c); //+= (char) c;

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
                for (LET i = i_lower; i < i_upper; ++i)
                    for (LET j = std::max(j_lower, i + 1); j < j_upper; ++j)
                        if (hall_set[j].first <= i) {
                            PARENT parents(i, j);
                            hall_set.push_back(parents);
                            degrees.push_back(d);
                            assert(((degrees[i] + degrees[j]) == degrees[hall_set.size() - 1]));
                            reverse_map[parents] = hall_set.size() - 1;
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
        assert (k <= size());
        return degrees[k];
    }

    /// Returns the key corresponding to a letter.
    inline KEY keyofletter(LET letter) const
    {
        return ltk.find(letter)->second;
    }

    /// Returns the left parent of a key.
    inline KEY lparent(const KEY &k) const
    {
        return hall_set[k].first;
    }

    /// Returns the right parent of a key.
    inline KEY rparent(const KEY &k) const
    {
        return hall_set[k].second;
    }

    /// Tells if a key corresponds to a letter.
    inline bool letter(const KEY &k) const
    {
        return ((k > 0) && (k <= letters.size()));
    }

    /// Returns the letter of a key corresponding to a letter.
    inline LET getletter(const KEY &k) const
    {
        return letters[k - 1];
    }

    /// Returns the value of the smallest key in the basis.
    inline KEY begin(void) const
    {
        return 1;
    }

    /// Returns the key next the biggest key of the basis.
    inline KEY end(void) const
    {
        return 0;
    }

    /// Returns the key next a given key in the basis. No implicit growup made.
    inline KEY nextkey(const KEY &k) const
    {
        if (k < (hall_set.size() - 1))
            return (k + 1);
        else
            return 0;
    }

    /// Returns the position of a key in the basis total order.
    inline DEG keypos(const KEY &k) const
    {
        return k;
    }

    /// Returns the size of the basis.
    inline KEY size(void) const
    {
        return (hall_set.size() - 1);
    }

    /// Outputs the Hall basis to an std::ostream.
    inline friend std::ostream &operator<<(std::ostream &os, const hall_basis &b)
    {
        for (KEY k = b.begin(); k != b.end(); k = b.nextkey(k))
            os << b.key2string(k) << ' ';
        return os;
    }

    //inline const std::string& key2string(const KEY& k) const
    //BUG//TJL//24/08/2012 - returned reference invalidated if vector grows!!
    //BUG//TJL//25/08/2012 - not templated but has static member so this is shared across all dimensions regardless of no letters etc!!
private:
    mutable std::vector<std::string> table; //move to instance per class
public:

    //ShortFix return a value not a reference
    //TODO check performance of fix 24/08/2012
    /// Converts a key to an std::string of letters.

    inline const std::string key2string(const KEY &k) const
    {
        static boost::recursive_mutex table_access;
        //// get exclusive recursive access for the thread
        boost::lock_guard<boost::recursive_mutex> lock(table_access);

        //BUG//TJL//09/04/2017 - non static member added to class but not commented out here!!
//		static std::vector<std::string> table;

        if (k > table.size()) {
            for (KEY i = (KEY) table.size() + 1; i <= k; ++i)
                table.push_back(_key2string(i));
        }
        return table[k - 1];
    }

private:
    /// Recursively constructs the string associated to the Lie key k.
    std::string _key2string(const KEY &k) const
    {
        std::ostringstream oss;
        if (k > 0) {
            if (letter(k))
                oss << getletter(k);
            else {
                oss << '[';
                oss << key2string(lparent(k));
                oss << ',';
                oss << key2string(rparent(k));
                oss << ']';
            }
        }
        return oss.str();
    }
};
//// if degree is static
//template<DEG n_letters>
//typename hall_basis<n_letters>::DEGREE hall_basis<n_letters>::degrees;

/// The Lie basis class.
/** 
 This is the basis used to implement the lie class as a specialisation of
 the algebra class. In the current implementation, the Lie basis class is a
 wrapper for the hall_basis class, with a prod() member function.
*/

template<typename SCA, typename RAT, DEG n_letters, DEG max_degree>
class lie_basis : protected hall_basis<n_letters>,
                  public basis_traits<With_Degree, n_letters, max_degree>,
                  dtl::hall_set_info<n_letters, max_degree>
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
    //using hall_basis<n_letters>::degrees;
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

    /// The MAP type.
    typedef std::map<KEY, SCA> MAP;
    /// The Free Lie Associative Algebra element type.
    struct temp_field
    {
        typedef SCA S;
        typedef RAT Q;
    };
    typedef lie<SCA, RAT, n_letters, max_degree,
        vectors::sparse_vector<lie_basis, temp_field> > LIE;
    /// The rationals.
    typedef SCA SCALAR;
    typedef RAT RATIONAL;

public:

    typedef basis::with_degree<max_degree> degree_tag;
    typedef basis::ordered<std::less<KEY> > ordering_tag;

public:
    /// Constructs the basis for a finite number of letters.
    lie_basis(void)
            : hall_basis<n_letters>()
    {
        // bug: tjl : 08 04 2017 without the following line the basis would not remain const and sharing it between threads would cause errors
        hall_basis<n_letters>::growup(max_degree);
    }

private:
    static std::map<PARENT, LIE> prime_prod_cache_table() {
        std::map<PARENT, LIE> rv;
        rv[PARENT(0, 0)] = LIE();
        return rv;
    }

public:
    /// Returns the product of two key.
    /**
    Returns the LIE instance corresponding to the product of the two basis
    elements k1 and k2. For performance reasons, the basis is enriched/grown
    dynamically and the already computed products are stored in a static
    multiplication table to speed up further calculations. This function
    returns a constant reference to the suitable table element. See Ch4 p93 and 94 Reutenauer
    */
    inline const LIE &prod(const KEY &k1, const KEY &k2)
    {
        static const LIE zero;
        DEG target_degree = degree(k1) + degree(k2); //degrees[k1] + degrees[k2];
        if ((max_degree > 0) && (target_degree > max_degree))
            return zero; // degree truncation
        // We grow up the basis up to the desired degree.
        growup(target_degree);

        static boost::recursive_mutex table_access;
        static std::map<PARENT, LIE> table(prime_prod_cache_table());
        // get exclusive recursive access for the thread
        boost::lock_guard<boost::recursive_mutex> lock(table_access);
        // [A,A] = 0.
        if (k1 == k2)
            return table[PARENT(0, 0)];
        typename std::map<PARENT, LIE>::iterator it;
        PARENT p(k1, k2);
        it = table.find(p);
        if (it == table.end()) {
            LIE *ptr =
                    &(table[p] =
                              ((p.first < p.second) ? _prod(k1, k2) : -_prod(k2, k1)));
            return *ptr;
        } else
            return it->second;
    }
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
        if (it != table.end())
            return it->second;
        else {
            if (letter(k)) {
                typename std::vector<LET>::size_type i;
                for (i = 0; i < s.size(); ++i)
                    if (s[i] == getletter(k))
                        return table[k] = *(v[i]);
                return (table[k] = (Lie) k);
            } else
                return (table[k]
                                = replace(lparent(k), s, v, table)
                                  * replace(rparent(k), s, v, table));
        }
    }

private:
    /// The recursive key product.
    // Caution Reutenauer only states the induction where the state variable
    // is the whole expression (the tree or the pre lie sum of products...)
    // we do not update that state and start over each time.
    LIE _prod(const KEY &k1, const KEY &k2)
    {
#ifdef DEBUG // NO LOOPS IN RECURSION FOR _PROD
        static std::map<PARENT, unsigned> counter;
        PARENT tmp(k1, k2);
        if (++counter[tmp] > 1) { assert(false); }; // add debug code here};
        assert(k1 < k2);
#endif // DEBUG

        // We look up for the desired product in our basis.
        PARENT parents(k1, k2);
        typename std::map<PARENT, KEY>::const_iterator it;
        it = reverse_map.find(parents);
        if (it != reverse_map.end())
            // [k1,k2] exists in the basis.
            return LIE(it->second);
        else
            // [k1,k2] does not exists in the basis.
        {
            // Since k1 <= k2, k2 is not a letter because if it was a letter,
            // then also k1, which is impossible since [k1,k2] is not in the basis.
            // similarly setting k2 = [k3,k4], it cannot be the case that k3 < k1
            // or this would be a hall tree. So k1 < k3 < k4
            // We use Jacobi: [k1,k2] = [k1,[k3,k4]]] = [[k1,k3],k4]-[[k1,k4],k3]
            KEY k3(lparent(k2));
            KEY k4(rparent(k2));
            LIE result(prod(k1, k3) * (LIE) k4);
            result.sub_mul(prod(k1, k4), (LIE) k3);
            return result;
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
    static DIMN key_to_index(const KEY k)
    {
        return DIMN(k - 1);
    }

    /// Convert an index to key
    static KEY index_to_key(const DIMN idx)
    {
        return KEY(idx + 1);
    }

    static DIMN start_of_degree(const DEG d)
    {
        assert(d <= max_degree + 1);
        return SIZE_INFO::degree_sizes[d];
    }


};


// Include once wrapper
#endif // DJC_COROPA_LIBALGEBRA_LIEBASISH_SEEN

//EOF.
