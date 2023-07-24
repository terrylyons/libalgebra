/* *************************************************************

Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurkó and Arend Janssen.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */

//  tensor_basis.h

// Include once wrapper
#ifndef DJC_COROPA_LIBALGEBRA_TENSORBASISH_SEEN
#define DJC_COROPA_LIBALGEBRA_TENSORBASISH_SEEN

#include "implementation_types.h"

#include <cmath>
#include <limits>

#include "_tensor_basis.h"
#include "base_basis.h"
#include "basis.h"
#include "detail/meta.h"
#include "key_iterators.h"

namespace alg {

namespace dtl {

using alg::integer_maths::power;

/*
 *  The second template argument is not really necessary, but we include this so
 *  that specializations can be partial rather than explicit. This means we can
 *  define statics inside the header rather than needing a cpp file.
 */

template<DEG NoLetters, typename I=DIMN>
struct tensor_size_info {
    static constexpr DEG bits_per_letter = ConstLog2<NoLetters - 1>::ans + 1;
    static constexpr DEG mantissa_bits_stored = std::numeric_limits<word_t>::digits - 1;
    static constexpr DEG max_depth = mantissa_bits_stored / bits_per_letter;

    template<std::size_t Depth>
    struct helper {
        static constexpr std::size_t value = integer_maths::sum_powers(NoLetters, Depth);
    };

    template<std::size_t D>
    struct power_helper {
        static constexpr std::size_t value = power(NoLetters, D);
    };

    using holder = typename alg::utils::generate_array<max_depth + 1, helper>::result;
    using power_holder = typename alg::utils::generate_array<max_depth, power_helper>::result;
    using degree_sizes_t = std::array<I, max_depth + 2>;

    static const std::array<I, max_depth + 1> powers;
    static const std::array<I, max_depth + 2> degree_sizes;
};

template<DEG N, typename I>
const std::array<I, tensor_size_info<N, I>::max_depth+ 2>
        tensor_size_info<N, I>::degree_sizes = tensor_size_info<N, I>::holder::data;

template<DEG Width, typename I>
const std::array<I, tensor_size_info<Width, I>::max_depth + 1> tensor_size_info<Width, I>::powers = tensor_size_info<Width, I>::power_holder::data;

template <typename I>
struct tensor_size_info<1, I> {
    static constexpr DEG max_depth = 52;

    template <std::size_t D>
    struct helper {
        static constexpr std::size_t value = D + 1;
    };

    template <std::size_t D>
    struct power_helper {
        static constexpr std::size_t value = 1;
    };

    using holder = typename alg::utils::generate_array<max_depth + 1, helper>::result;
    using power_holder = typename alg::utils::generate_array<max_depth, power_helper>::result;
    using degree_sizes_t = std::array<I, max_depth + 2>;

    static const std::array<I, max_depth+1> powers;
    static const std::array<I, max_depth+2> degree_sizes;
};


template <typename I>
const std::array<I, tensor_size_info<1, I>::max_depth + 2>
        tensor_size_info<1, I>::degree_sizes = tensor_size_info<1, I>::holder::data;

template <typename I>
const std::array<I, tensor_size_info<1, I>::max_depth + 1> tensor_size_info<1, I>::powers = tensor_size_info<1, I>::power_holder::data;


}// namespace dtl

/// Implements an interface for the set of words of a finite number of letters.
/**

A basis is a finite total ordered set of keys, its cardinal is size() and
its minimal element is begin(). The successor key of a given key is given
by nextkey(). The successor of the maximal key is end() and does not
belong to the basis. The position of a given key in the total order of the
basis is given by keypos(), and equals 1 for begin(). To each letter
corresponds a key.

The tensor_basis is a basis for which keys are words of letters, totally
ordered by the lexicographical order. Letters corresponds to words of
length one. The empty word is a special key which serves for the embedding
of scalars into the tensor algebra. The begin() key is the empty word.

Implementation: a key is implemeted using dual ended queue of letters.
Letters are natural numbers, enumerated starting from 1, in such a way
that n_letters is the biggest letter. We use the lexicographical order.
The invalid key made with one occurrence of letter 0 is used for end().

There is no key product here since the tensor_basis class serves as a
common ancestor for free_tensor_basis and shuffle_tensor_basis classes.
The implementation of the prod() member function constitutes the essential
difference between free_tensor_basis and shuffle_tensor_basis.
*/
template<DEG n_letters, DEG max_degree>
class tensor_basis : dtl::tensor_size_info<n_letters>
{
public:

    typedef dtl::tensor_size_info<n_letters> SIZE_INFO;
    /// A key is an dual ended queue of letters
    typedef _tensor_basis<n_letters, max_degree> KEY;
    /// A default key corresponds to the empty word.
    const KEY empty_key;
    /// The MAP type.
#if 0
#ifndef ORDEREDMAP
    typedef MY_UNORDERED_MAP <KEY, SCA, typename KEY::hash> MAP;
#else
#ifdef NOBTREE
    typedef std::map<KEY, SCA> MAP;
#else
    typedef btree::safe_btree_map<KEY, SCA, std::less<KEY>, std::allocator<std::pair<const KEY, SCA> >, 256> MAP;
#endif
#endif// !ORDEREDMAP
#endif
    typedef alg::basis::with_degree<max_degree> degree_tag;
    typedef alg::basis::ordered<std::less<KEY>> ordering_tag;

    static const DEG s_no_letters = n_letters;
    static const DEG s_max_degree = max_degree;

private:
    /// The size of the full basis
    LET _size;

public:
    /// Default constructor. Empty basis.

    /// We statically compute the size (n^{d+1}-1)/(n-1)
    /// where n == n_letters and d == max_degree.
    tensor_basis(void)
        : _size(dtl::tensor_size_info<n_letters>::degree_sizes[max_degree])
    {}

public:
    /// Returns the key corresponding to a letter.
    inline KEY keyofletter(LET letter) const
    {
        return KEY(letter);
    }

    /// Tells if a key is a letter (i.e. word of length one).
    inline bool letter(const KEY& k) const
    {
        return k.size() == 1;
    }

    /// Returns the first letter of a key.
    inline LET getletter(const KEY& k) const
    {
        return k.FirstLetter();
    }

    /// Returns the first letter of a key. For compatibility with lie_basis.
    inline KEY lparent(const KEY& k) const
    {
        return k.lparent();
    }

    /// Returns the key which corresponds to the sub-word after the first letter.
    inline KEY rparent(const KEY& k) const
    {
        return k.rparent();
    }

    /// Returns the length of the key viewed as a word of letters.
    inline DEG degree(const KEY& k) const
    {
        return k.size();
    }

    /// Returns the size of the basis.
    inline LET size(void) const
    {
        return LET(_size);
    }

    /// Returns the value of the smallest key in the basis.
    inline static KEY begin(void)
    {
        return KEY();
    }

    /// Returns the key next the biggest key of the basis.
    inline static KEY end(void)
    {
        // TJL 21/08/2012
        // KEY result; // empty key.
        // result.push_back(0); // invalid key.
        return KEY::end();
    }

    /// Returns the key next a given key in the basis.
    inline static KEY nextkey(const KEY& k)
    {
        // tjl  25/08/2012
        size_t i = k.size();
        KEY result(k);
        for (size_t j = 0; j < i; ++j) {
            if (k[i - 1 - j] < LET(n_letters)) {
                result[i - 1 - j] += 1;
                return result;
            }
            else {
                result[i - 1 - j] = KEY(LET(1)).FirstLetter();
            }
        }
        if (k.size() == max_degree) {
            return end();
        }
        else {
            return KEY(LET(1)) * result;
        }
    }

    /// Outputs a key as a string of letters to an std::ostringstream.
    std::string key2string(const KEY& k) const
    {
        std::ostringstream oss;

        if (k.size() > 0) {
            oss << k.FirstLetter();
            KEY kk = k.rparent();
            for (unsigned i = 1; i < k.size(); ++i) {
                oss << "," << kk.FirstLetter();
                kk = kk.rparent();
            }
        }
        return oss.str();
    }

    // Key iteration methods

    basis::key_range<tensor_basis> iterate_keys() const noexcept
    {
        return basis::key_range<tensor_basis>(*this);
    }

    basis::key_range<tensor_basis> iterate_keys_to_deg(DEG mdeg) const noexcept
    {
        if (mdeg > max_degree) {
            return iterate_keys();
        }
        return basis::key_range<tensor_basis>(*this, KEY(), index_to_key(start_of_degree(mdeg)));
    }

    basis::key_range<tensor_basis> iterate_keys(const KEY& begin, const KEY& end) const noexcept
    {
        return basis::key_range<tensor_basis>{*this, begin, end};
    }

    basis::key_range<tensor_basis> iterate_keys_from(const KEY& begin) const noexcept
    {
        return basis::key_range<tensor_basis>{*this, begin};
    }

    basis::key_range<tensor_basis> iterate_keys_to(const KEY& end) const noexcept
    {
        return basis::key_range<tensor_basis>{*this, begin(), end};
    }

private:
    /*
    static DIMN key_to_index_impl(KEY const& key) {

    }
*/
public:
    /*    static DIMN key_to_index(const KEY &key)
    {
        //assert(key.valid());
        DEG size = key.size();

        if (size == 0) {
            return 0;
        } else if (size == 1) {
            return static_cast<DIMN>(key.FirstLetter());
        }

        static std::map<KEY, DIMN> table = []() {
            std::map<KEY, DIMN> t;
            KEY k(LET(1));
            k.push_back(LET(1));

            DIMN idx = n_letters;
            while (k.size() == 2) {
                t[k] = ++idx;
                k = nextkey(k);
            }
            return t;
        }();


        DEG s2 = size / 2;
        KEY right(key), left=right.split_n(s2);
        assert(left.size() == s2);
        assert(right.size() == (size - s2));
        return key_to_index(left) * pow(n_letters, (size-s2)) + key_to_index(right);


    }
*/

    static DIMN key_to_index(const KEY& key)
    {
        assert(key.valid());

        DIMN idx = 0;
        if (key.size() == 0) {
            return idx;
        }

        KEY tmp(key);
        while (tmp.size() > 0) {
            idx *= n_letters;
            idx += static_cast<DIMN>(tmp.FirstLetter());
            tmp = tmp.rparent();
        }
        return idx;
    }

    static KEY index_to_key(const DIMN idx)
    {
        if (idx == 0) {
            return KEY();
        }
        else if (1 <= idx && idx <= n_letters) {
            return KEY(LET(idx));
        }

//        static boost::recursive_mutex access;
//        boost::lock_guard<boost::recursive_mutex> lock(access);
        // static const std::vector<KEY> __cache = _key_of_index_cache();
        // assert(idx < __cache.size());
        // return __cache[idx];

        //static std::map<DIMN, KEY> cache;

        //typename std::map<DIMN, KEY>::iterator it = cache.find(idx);

        //if (it != cache.end()) {
        //    return it->second;
        //}

//        std::map<DIMN, KEY> cache;
//
//        KEY& rv = cache[idx];
//
//        if (rv.size() > 0) {
//            return rv;
//        }

        DIMN i = idx;

        DIMN l;
        KEY val;
        while (i) {
            i -= 1;
            l = i % n_letters;
            val.push_back(LET(1 + l));
            i /= n_letters;
        }

        return /*rv =*/ val.reverse();
    }

    static DIMN start_of_degree(const DEG deg)
    {
        assert(deg <= max_degree + 1);
        return (deg == 0) ? 0 : SIZE_INFO::degree_sizes[deg - 1];
    }

public:
    /// Outputs a std::pair<shuffle_tensor_basis*, KEY> to an std::ostream.
    inline friend std::ostream& operator<<(std::ostream& os, const std::pair<tensor_basis*, KEY>& t)
    {
        return os << (t.first)->key2string(t.second);
    }
};

/**
 * @brief The monoid of words of a finite number of letters with concat product.
 *
 * This is the basis used to implement the free_tensor class as a
 * specialisation of the algebra class. This basis is the Free Associative
 * Algebra basis with a finite number of letters, with the usual
 * concatenation product. The free_tensor_basis is a container of keys. A key
 * is the implementation of a word of letters. The prod() member function
 * corresponds to the concatenation of the two keys given as arguments. This
 * product is associative but not commutative. Letters can be seen as
 * particular basis keys, i.e. words of length one. The empty word is a
 * special key used for the imbedding of letters (words of length one).
*/
template<DEG n_letters, DEG max_degree>
class free_tensor_basis : public tensor_basis<n_letters, max_degree>,
                          public base_basis<With_Degree, n_letters, max_degree>
{
public:
    /// The tensor_basis type.
    typedef tensor_basis<n_letters, max_degree> TBASIS;

    typedef typename TBASIS::KEY KEY;

    /// The Free Associative Algebra element type.

public:
    /// Default constructor.
    free_tensor_basis(void)
    {}

public:
    /// The concatenation product of two basis elements.
    /**
    Returns the free_tensor obtained by the concatenation product of two keys
    viewed as words of letters. The result is a unidimensional free_tensor
    with a unique key (the concatenation of k1 and k2) associated to the +1
    scalar. The already computed products are not stored or remembered.
    */

    // static inline TENSOR prod(const typename alg::tensor_basis<SCA, n_letters,
    // max_degree>::KEY& k1, 	const typename alg::tensor_basis<SCA, n_letters,
    //max_degree>::KEY& k2)
    //{
    //	SCA one(+1);
    //	TENSOR result;
    //	if ((max_degree == 0) || (k1.size() + k2.size() <= max_degree))
    //	{
    //		//typename alg::tensor_basis<SCA, n_letters, max_degree>::KEY concat
    //= k1 * k2; 		result[k1 * k2] = one;
    //	}
    //	else
    //	{
    //		throw;
    //	}
    //	return result;
    // }

    static inline typename TBASIS::KEY prod(const typename TBASIS::KEY& k1, const typename TBASIS::KEY& k2)
    {
        assert((max_degree == 0) || (k1.size() + k2.size() <= max_degree));
        return k1 * k2;
    }

    // static inline TENSOR& prod(const typename alg::tensor_basis<SCA, n_letters,
    // max_degree>::KEY& k1, 	const typename alg::tensor_basis<SCA, n_letters,
    //max_degree>::KEY& k2, 	TENSOR& ans = (TENSOR()))
    //{
    //	assert((max_degree == 0) || (k1.size() + k2.size() <= max_degree));
    //	ans[k1 * k2] = SCA(+1);
    //	return ans;
    // }

    /// Outputs an std::pair<free_tensor_basis*, KEY> to an std::ostream.
    inline friend std::ostream& operator<<(std::ostream& os, const std::pair<free_tensor_basis*, KEY&>& t)
    {
        return os << (t.first)->key2string(t.second);
    }
};

namespace vectors {

template<DEG n_letters, DEG max_depth, typename Field>
struct vector_type_selector<free_tensor_basis<n_letters, max_depth>, Field> {
    typedef free_tensor_basis<n_letters, max_depth> BASIS;
    typedef typename BASIS::KEY KEY;
    typedef vectors::sparse_vector<BASIS, Field,
#ifndef ORDEREDMAP
                                   MY_UNORDERED_MAP<KEY, typename Field::S, typename KEY::hash>
#else
                                   std::map<KEY, typename Field::S>
#endif
                                   >
            sparse_vect;

    typedef vectors::dense_vector<BASIS, Field> dense_vect;

    typedef vectors::hybrid_vector<BASIS, Field, vectors::policy::basic_resize_policy,
#ifndef ORDEREDMAP
                                   MY_UNORDERED_MAP<KEY, typename Field::S, typename KEY::hash>
#else
                                   std::map<KEY, typename Field::S>
#endif
                                   >
            hybrid_vect;

    typedef typename alg::utils::type_selector<boost::is_pod<typename Field::S>::value, sparse_vect, dense_vect>::type
            type;
};

}// namespace vectors

/**
 * @brief The monoid of words of a finite number of letters with shuffle product.
 *
 * This is the basis used to implement the shuffle_tensor class as a
 * specialisation of the algebra class. This basis is the Free Associative
 * Algebra basis with a finite number of letters, with the shuffle product.
 * The shuffle_tensor_basis is a container of keys. A key is the
 * implementation of a word of letters. The prod() member function
 * corresponds to the shuffle product of the two keys given as arguments.
 * This product is associative and commutative. Letters can be seen as
 * particular basis keys, i.e. words of length one. The empty word is a
 * special key used for the embedding of letters (words of length one).
*/
template<DEG n_letters, DEG max_degree>
class shuffle_tensor_basis : public tensor_basis<n_letters, max_degree>,
                             public base_basis<With_Degree, n_letters, max_degree>
{
public:
    /// The tensor_basis type.
    typedef tensor_basis<n_letters, max_degree> TBASIS;
    /// Import of the KEY type.
    typedef typename TBASIS::KEY KEY;
    /// Import of the MAP type.
    // typedef typename TBASIS::MAP MAP;

    typedef alg::basis::with_degree<max_degree> degree_tag;
    typedef alg::basis::ordered<std::less<KEY>> ordering_tag;

public:
    /// Default constructor.
    shuffle_tensor_basis(void)
    {}
};
//#endif

namespace vectors {

template<DEG n_letters, DEG max_depth, typename Field>
struct vector_type_selector<shuffle_tensor_basis<n_letters, max_depth>, Field> {
    typedef shuffle_tensor_basis<n_letters, max_depth> BASIS;
    typedef typename BASIS::KEY KEY;
    typedef sparse_vector<BASIS, Field,
#ifndef ORDEREDMAP
                          MY_UNORDERED_MAP<KEY, typename Field::S, typename KEY::hash>
#else
                          std::map<KEY, typename Field::S>
#endif
                          >
            type;
};

}// namespace vectors

namespace basis {

template<DEG Width, DEG Depth1, DEG Depth2>
struct related_to<tensor_basis<Width, Depth1>, tensor_basis<Width, Depth2>>
    : std::true_type {
};

template<DEG Width, DEG Depth1, DEG Depth2>
struct related_to<free_tensor_basis<Width, Depth1>, free_tensor_basis<Width, Depth2>>
    : std::true_type {
};

template<DEG Width, DEG Depth1, DEG Depth2>
struct related_to<shuffle_tensor_basis<Width, Depth1>, shuffle_tensor_basis<Width, Depth2>>
    : std::true_type {
};

template<DEG Width, DEG Depth1, DEG Depth2>
struct related_to<free_tensor_basis<Width, Depth1>, tensor_basis<Width, Depth2>>
    : std::true_type {
};

template<DEG Width, DEG Depth1, DEG Depth2>
struct related_to<tensor_basis<Width, Depth1>, free_tensor_basis<Width, Depth2>>
    : std::true_type {
};

template<DEG Width, DEG Depth1, DEG Depth2>
struct related_to<shuffle_tensor_basis<Width, Depth1>, tensor_basis<Width, Depth2>>
    : std::true_type {
};

template<DEG Width, DEG Depth1, DEG Depth2>
struct related_to<tensor_basis<Width, Depth1>, shuffle_tensor_basis<Width, Depth2>>
    : std::true_type {
};

}// namespace basis

}// namespace alg
// Include once wrapper
#endif// DJC_COROPA_LIBALGEBRA_TENSORBASISH_SEEN

// EOF.
