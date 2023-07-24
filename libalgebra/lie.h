/* *************************************************************

Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurkó and Arend Janssen.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */

//  lie.h

// Include once wrapper
#ifndef DJC_COROPA_LIBALGEBRA_LIEH_SEEN
#define DJC_COROPA_LIBALGEBRA_LIEH_SEEN

#include "algebra.h"
#include "lie_basis.h"



namespace alg {


template <DEG Width, DEG Depth>
class lie_multiplier
    : public multiplier_base<lie_multiplier<Width, Depth>,
        lie_basis<Width, Depth>, 2>
{
    using base = multiplier_base<lie_multiplier<Width, Depth>, lie_basis<Width, Depth>, 2>;
    friend base;

    static const lie_basis<Width, Depth> basis;

public:

    using basis_type = lie_basis<Width, Depth>;
    using key_type = typename basis_type::KEY;
    using scalar_type = int;
    using typename base::pair_type;
    using typename base::inner_result_type;
    using typename base::result_type;
    using typename base::argument_type;

private:

    inner_result_type prod_impl(argument_type lhs, argument_type rhs) const
    {
        if (rhs < lhs) {
            return base::uminus(prod_impl(rhs, lhs));
        }

        auto found = basis.find({lhs, rhs});
        if (found != basis.reverse_map_end()) {
            return {{found->second, 1}};
        }

        auto lparent = basis.lparent(rhs);
        auto rparent = basis.rparent(rhs);

        auto left_result = base::mul(operator()(lhs, lparent), rparent);
        auto right_result = base::mul(operator()(lhs, rparent), lparent);

        return base::sub(left_result, right_result);
    }

public:

    result_type operator()(argument_type lhs, argument_type rhs) const
    {
        static const boost::container::small_vector<pair_type, 0> null;

        if (basis.degree(lhs) + basis.degree(rhs) > Depth) {
            return null;
        }
        if (lhs == rhs) {
            return null;
        }

        return base::cached_compute(lhs, rhs);
    }



};

template <DEG Width, DEG Depth>
const lie_basis<Width, Depth> lie_multiplier<Width, Depth>::basis;

template <DEG Width, DEG Depth>
using lie_multiplication = base_multiplication<lie_multiplier<Width, Depth>>;


/**
 * @brief  A specialisation of the algebra class with a Lie basis.
 *
 * Mathematically, the algebra of Lie instances is a free Lie associative
 * algebra. With respect to the inherited algebra class, the essential
 * distinguishing feature of this class is the basis class used, and in
 * particular the basis::prod() member function. Thus, the most important
 * information is in the definition of lie_basis. Notice that this associative
 * algebra of lie elements does not includes as a sub-algebra the associative
 * algebra corresponding to the SCALAR type. In other words, only the scalar
 * zero corresponds to a Lie element (the zero one) which is the neutral
 * element of the addition operation. There is no neutral element for the
 * product (free Lie product).
 */
template<typename Coeff, DEG n_letters, DEG max_degree,
         template<typename, typename, typename...> class VectorType,
         typename... Args>
class lie : public algebra<
                    lie_basis<n_letters, max_degree>, Coeff, lie_multiplication<n_letters, max_degree>, VectorType, Args...>
{
    typedef lie_multiplication<n_letters, max_degree> multiplication_t;

public:
    /// The basis type.
    typedef lie_basis<n_letters, max_degree> BASIS;
    /// Import of the KEY type.
    typedef typename BASIS::KEY KEY;
    /// The algebra type.
    typedef algebra<BASIS, Coeff, multiplication_t, VectorType, Args...> ALG;
    /// The sparse_vector type.
    typedef typename ALG::VECT VECT;

    typedef typename Coeff::SCA SCA;
    typedef typename Coeff::RAT RAT;

    /// Import of the iterator type.
    typedef typename ALG::iterator iterator;
    /// Import of the constant iterator type.
    typedef typename ALG::const_iterator const_iterator;

public:
    /// Default constructor. Zero lie element.
    lie()
    {}

    /// Copy constructor.
    lie(const lie& l)
        : ALG(l)
    {}

    /// Constructs an instance from an algebra instance.
    lie(const ALG& a)
        : ALG(a)
    {}

    /// Constructs an instance from a sparse_vector instance.
    lie(const VECT& v)
        : ALG(v)
    {}

    /// Constructs a unidimensional instance from a given key (with scalar one).
    explicit lie(const KEY& k)
        : ALG(k)
    {}

    // explicit lie(LET letter, const SCA& s)
    //	// flawed as basis is possibly not yet constructed
    //	: ALG(VECT::basis.keyofletter(letter), s) {}
    /// Constructs a unidimensional instance from a key and a scalar.
    explicit lie(const KEY& k, const SCA& s)
        : ALG(k, s)
    {}

    /// Construct an instance from pointers to data range
    lie(SCA const* begin,
        SCA const* end) : ALG(begin, end)
    {
    }

    template <typename InputIt>
    lie(InputIt begin, InputIt end) : ALG(begin, end)
    {}

    lie& operator=(const lie&) = default;
    //lie& operator=(lie&&) noexcept = default;

public:
    /// Replaces the occurrences of letters in s by Lie elements in v.
    inline friend lie replace(const lie& src, const std::vector<LET>& s, const std::vector<const lie*>& v)
    {
        lie result;
        std::map<KEY, lie> table;
        const_iterator i;
        for (i = src.begin(); i != src.end(); ++i) {
            result.add_scal_prod(VECT::basis.replace(i->key(), s, v, table), i->value());
        }
        return result;
    }


#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
private:

    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, unsigned int const /* version */) {
        ar & boost::serialization::base_object<ALG>(*this);
    }
#endif
};

}// namespace alg

// Include once wrapper
#endif// DJC_COROPA_LIBALGEBRA_LIEH_SEEN

// EOF.
