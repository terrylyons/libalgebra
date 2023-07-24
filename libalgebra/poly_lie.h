/* *************************************************************

Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurkó and Arend Janssen.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */

#ifndef DJC_COROPA_LIBALGEBRA_POLYLIEH_SEEN
#define DJC_COROPA_LIBALGEBRA_POLYLIEH_SEEN

#include "poly_basis.h"
#include "poly_lie_basis.h"
#include "polynomials.h"


namespace alg {


template <DEG Width, DEG Depth>
class poly_lie_multiplier : public multiplier_base<poly_lie_multiplier<Width, Depth>,
        poly_lie_basis<Width, Depth>, 2>
{
    using base = multiplier_base<poly_lie_multiplier<Width, Depth>,
            poly_lie_basis<Width, Depth>, 2>;
    friend base;

    poly_multiplier polymul;

public:

    using basis_type = poly_lie_basis<Width, Depth>;
    using key_type = typename basis_type::KEY;
    using result_type = typename base::inner_result_type;
    using typename base::argument_type;

private:
    using poly_result_type = typename poly_multiplier::result_type;
    using poly_reference = const boost::container::small_vector_base<
            typename poly_multiplier::pair_type>&;
    using poly_argument_type = typename poly_multiplier::argument_type;
    using poly_let_type = typename poly_basis::LET;
    using poly_key_type = typename poly_basis::KEY;

    result_type prod2(poly_result_type pol, argument_type key) const
    {
        result_type result;
        result.reserve(pol.size());
        for (const auto& item : pol) {
            result.emplace_back(
                    key_type(key.first, polymul(key.second, item.first)[0].first),
                    item.second
                    );
        }
        return result;
    }

    poly_result_type partial_diff(poly_argument_type mon, poly_let_type letter) const
    {
        poly_result_type result;

        poly_key_type mon_c(mon);
        auto found = mon_c.find(letter);
        if (found != mon_c.end()) {
            if (found->second == 1) {
                mon_c.erase(found);
                result.emplace_back(mon_c, 1);
            } else {
                auto coeff = (found->second)--;
                result.emplace_back(mon_c, coeff);
            }
        }
        return result;
    }

public:

    result_type operator()(argument_type lhs, argument_type rhs) const
    {
        auto poly1 = partial_diff(rhs.second, lhs.first);
        auto poly2 = partial_diff(lhs.second, rhs.first);
        key_type mon1(rhs.first, lhs.second);
        key_type mon2(lhs.first, rhs.second);
        return base::sub(prod2(poly1, mon1), prod2(poly2, mon2));
    }
};


template <DEG Width, DEG Depth>
using poly_lie_multiplication = base_multiplication<poly_lie_multiplier<Width, Depth>>;


/// The Lie algebra for the commutative polynomials.

/// Elements of the algebra are polynomial vector fields (ie linear combinations
/// of pairs of monomials and directions). The product is the Lie bracket
/// for vector fields.

template<typename Coeff, DEG n_letters, DEG max_degree, template<typename, typename, typename...> class VectorType, typename... Args>
class poly_lie : public algebra<
                         poly_lie_basis<n_letters, max_degree>,
                         Coeff,
                         poly_lie_multiplication<n_letters, max_degree>,
                         VectorType, Args...>
{
    typedef poly_lie_multiplication<n_letters, max_degree> multiplication_t;

public:
    typedef typename Coeff::S SCA;
    typedef typename Coeff::Q RAT;

    /// The basis elements for the polynomials (monomials)
    typedef typename poly_basis::KEY POLYBASISKEY;
    /// The basis type.
    typedef poly_lie_basis<n_letters, max_degree> BASIS;
    /// Import of the KEY type.
    typedef typename BASIS::KEY KEY;
    /// The algebra type.
    typedef algebra<BASIS, Coeff, multiplication_t, VectorType, Args...> ALG;
    /// The sparse_vector type.
    typedef typename ALG::VECT VECT;
    /// Import of the iterator type.
    typedef typename ALG::iterator iterator;
    /// Import of the constant iterator type.
    typedef typename ALG::const_iterator const_iterator;

public:
    /// Default constructor. Zero lie element.
    poly_lie()
    {}

    /// Copy constructor.
    poly_lie(const poly_lie& l)
        : ALG(l)
    {}

    /// Constructs an instance from an algebra instance.
    poly_lie(const ALG& a)
        : ALG(a)
    {}

    /// Constructs an instance from a sparse_vector instance.
    poly_lie(const VECT& v)
        : ALG(v)
    {}

    /// Constructs a unidimensional instance from a given key (with scalar one).
    explicit poly_lie(LET x, LET y, DEG z)
        : ALG()
    {
        POLYBASISKEY tempkey;
        tempkey[y] = z;
        ALG tempalg(KEY(x, tempkey));
        ALG::swap(tempalg);
    }

    /// Constructs an instance from a basis element.
    explicit poly_lie(const KEY& k)
        : ALG(k)
    {}

    /// Constructs a unidimensional instance from a letter and a scalar.
    explicit poly_lie(LET letter, const SCA& s)
        : ALG(VECT::basis.keyofletter(letter), s)
    {}

    poly_lie& operator=(const poly_lie&) = default;
    poly_lie& operator=(poly_lie&&) noexcept = default;

public:
    /// Replaces the occurrences of letters in s by Lie elements in v.
    inline friend poly_lie
    replace(const poly_lie& src, const std::vector<LET>& s, const std::vector<const poly_lie*>& v)
    {
        poly_lie result;
        std::map<KEY, poly_lie> table;
        const_iterator i;
        for (i = src.begin(); i != src.end(); ++i) {
            result.add_scal_prod(VECT::basis.replace(i->first, s, v, table), i->second);
        }
        return result;
    }

#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
private:
    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive& ar, unsigned int const /* version */)
    {
        ar& boost::serialization::base_object<ALG>(*this);
    }
#endif
};

}// namespace alg

#endif
