/* *************************************************************

Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurk� and Arend Janssen.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */

#ifndef DJC_COROPA_LIBALGEBRA_POLYLIEH_SEEN
#define DJC_COROPA_LIBALGEBRA_POLYLIEH_SEEN

#include "libalgebra/polynomials.h"

namespace alg {

template<typename Coeff>
class poly_lie_multiplication : poly_multiplication<Coeff>
{
    typedef typename Coeff::SCA scalar_t;
    typedef poly_multiplication<Coeff> poly_multiplication_t;

    /// Returns the Lie bracket of two monomial vector fields.
    /**
    Returns the Vector field corresponding to the Lie bracket of two vector
    fields. If we have have monomials m1 and m2 and Letters i and j then the
    product of m1*d/dxi with m2*d/dxj is equal to m1*(d/dxi m2)*d/dxj - m2*(d/dxj
    m1)*d/dxi
    */
    template<typename PolyLie>
    static PolyLie prod(typename PolyLie::KEY const& k1, typename PolyLie::KEY const& k2)
    {
        typedef poly<Coeff> poly_t;
        typedef typename PolyLie::KEY key_t;
        poly_t poly1 = poly_t::prediff(k2.second, k1.first);
        poly_t poly2 = poly_t::prediff(k1.second, k2.first);
        key_t mon1(k2.first, k1.second);
        key_t mon2(k1.first, k2.second);
        PolyLie result;
        result = prod2<PolyLie>(poly1, mon1) - prod2<PolyLie>(poly2, mon2);
        return result;
    }

    /// Multiplication of a polynomial poly1 by a monomial vector field liemon1.
    template<typename PolyLie>
    static PolyLie prod2(poly<Coeff> const& poly1, typename PolyLie::KEY const& liemon1)
    {
        typedef poly<Coeff> poly_t;
        //typedef typename PolyLie::KEY key_t;
        PolyLie result;
        for (typename poly_t::const_iterator it = poly1.begin(); it != poly1.end(); it++) {
            scalar_t temp = it->value();
            typename poly_t::BASIS::KEY temp2 = poly_multiplication_t::prod2(liemon1.second, it->key());
            result[make_pair(liemon1.first, temp2)] = temp;
        }
        return result;
    }

    template<typename Transform, typename Algebra>
    class key_operator
    {

        Transform m_transform;

    public:
#if __cplusplus >= 201103UL

        template<typename... Args>
        explicit key_operator(Args&&... args)
            : m_transform(std::forward<Args>(args)...)
        {}

#else
        template<typename Arg>
        explicit key_operator(Arg arg) : m_transform(arg)
        {}
#endif

        template<typename Vector>
        void
        operator()(Vector& result, typename Vector::KEY const& lhs_key, scalar_t const& lhs_val,
                   typename Vector::KEY const& rhs_key, scalar_t const& rhs_val)
        {
            result.add_scal_prod(prod<Algebra>(lhs_key, rhs_key), m_transform(Coeff::mul(lhs_val, rhs_val)));
        }
    };

public:
    template<typename Algebra, typename Operator>
    Algebra& multiply_and_add(Algebra& result, Algebra const& lhs, Algebra const& rhs, Operator op) const
    {
        key_operator<Operator, Algebra> kt(op);
        lhs.buffered_apply_binary_transform(result, rhs, kt);
        return result;
    }

    template<typename Algebra, typename Operator>
    Algebra&
    multiply_and_add(Algebra& result, Algebra const& lhs, Algebra const& rhs, Operator op, DEG const max_depth) const
    {
        key_operator<Operator, Algebra> kt(op);
        lhs.buffered_apply_binary_transform(result, rhs, kt, max_depth);
        return result;
    }

    template<typename Algebra, typename Operator>
    Algebra multiply(Algebra const& lhs, Algebra const& rhs, Operator op) const
    {
        Algebra result;
        multiply_and_add(result, lhs, rhs, op);
        return result;
    }

    template<typename Algebra, typename Operator>
    Algebra multiply(Algebra const& lhs, Algebra const& rhs, Operator op, DEG const max_depth) const
    {
        Algebra result;
        multiply_and_add(result, lhs, rhs, op, max_depth);
        return result;
    }

    template<typename Algebra, typename Operator>
    Algebra& multiply_inplace(Algebra& lhs, Algebra const& rhs, Operator op) const
    {
        key_operator<Operator, Algebra> kt(op);
        lhs.unbuffered_apply_binary_transform(rhs, kt);
        return lhs;
    }

    template<typename Algebra, typename Operator>
    Algebra& multiply_inplace(Algebra& lhs, Algebra const& rhs, Operator op, DEG const max_depth) const
    {
        key_operator<Operator, Algebra> kt(op);
        lhs.unbuffered_apply_binary_transform(rhs, kt, max_depth);
        return lhs;
    }
};

/// The Lie algebra for the commutative polynomials.

/// Elements of the algebra are polynomial vector fields (ie linear combinations
/// of pairs of monomials and directions). The product is the Lie bracket
/// for vector fields.

template<typename Coeff, DEG n_letters, DEG max_degree, typename...>
class poly_lie : public algebra<
                         poly_lie_basis<n_letters, max_degree>,
                         Coeff,
                         poly_lie_multiplication<Coeff>>
{
    typedef poly_lie_multiplication<Coeff> multiplication_t;

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
    typedef algebra<BASIS, Coeff, multiplication_t> ALG;
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
    replace(const poly_lie &src, const std::vector<LET> &s, const std::vector<const poly_lie *> &v)
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

template <typename Archive>
void serialize(Archive &ar, unsigned int const /* version */) {
    ar & boost::serialization::base_object<ALG>(*this);
}
#endif
};

}// namespace alg

#endif
