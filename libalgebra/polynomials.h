/* *************************************************************

Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurkó and Arend Janssen.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */

//  polynomials.h

// Include once wrapper
#ifndef DJC_COROPA_LIBALGEBRA_POLYNOMIALSH_SEEN
#define DJC_COROPA_LIBALGEBRA_POLYNOMIALSH_SEEN

#include "algebra.h"
#include "poly_basis.h"


namespace alg {


class poly_multiplier : public multiplier_base<poly_multiplier, poly_basis>
{
    using base = multiplier_base<poly_multiplier, poly_basis>;
    friend base;

public:
    using basis_type = poly_basis;
    using key_type = typename poly_basis::KEY;
    using typename base::pair_type;
    using result_type = typename base::inner_result_type;
    using typename base::argument_type;


    result_type operator()(argument_type lhs, argument_type rhs) const
    {
        key_type result(lhs);
        for (const auto& item : rhs) {
            result[item.first] += item.second;
        }
        return {{result, 1}};
    }

};


using poly_multiplication = base_multiplication<poly_multiplier>;


/// A specialisation of the algebra class with a commutative monomial product.
/**
   Mathematically, the algebra of polynomials instances is a communtative and
   associative polynomial algebra in several variables (letters). With respect
   to the inherited algebra class, the essential distinguishing feature of
   this class is the basis class used, and in particular the basis::prod()
   member function. Thus, the most important information is in the definition
   of poly_basis. Notice that this associative algebra of polynomials includes
   as a sub-algebra the associative algebra corresponding to the SCALAR type.
   This is permitted by the existence of empty keys in poly_basis. A
   polynomial is essentially a sparse vector of monomials with polynomial
   commutative product.
 */
template<typename Coeff, typename...>
class poly : public algebra<poly_basis, Coeff, poly_multiplication>
{
    typedef poly_multiplication multiplication_t;
public:
    typedef typename Coeff::S SCA;
    typedef typename Coeff::Q RAT;

    /// The basis type.
    typedef poly_basis BASIS;
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

    /// Default constructor. Empty polynomial. Zero.
    poly() = default;

    /// Copy constructor.
    poly(const poly& p)
        : ALG(p)
    {}

    /// Constructs an instance from an algebra instance.
    poly(const ALG& a)
        : ALG(a)
    {}

    /// Constructs an instance from a sparse_vector instance.
    poly(const VECT& v)
        : ALG(v)
    {}

    /// Constructs an instance from a scalar. Embedding of scalars.
    explicit poly(const SCA& s)
        : ALG(BASIS::KEY(), s)
    {}

    /// Constructs a unidimensional instance from a key (a monomial).
    explicit poly(const KEY& k)
        : ALG(k)
    {}

    /// Constructs a unidimensional instance from a given letter and scalar.
    explicit poly(LET letter, const SCA& s)
        : ALG(VECT::basis.keyofletter(letter), s)
    {}

    poly& operator=(const poly&) = default;
    poly& operator=(poly&&) noexcept = default;

    /// Evaluates the polynomial for some scalar values for letters
    /// (variables).
    inline SCA eval(const std::map<LET, SCA>& values) const
    {
        SCA result(VECT::zero);
        for (const_iterator i = VECT::begin(); i != VECT::end(); ++i) {
            Coeff::add_inplace(result, Coeff::mul(VECT::basis.template eval_key<Coeff>(i->key(), values), i->value()));
        }
        return result;
    }

public:
    /// Partial differentiation of the KEY (monomial) k1 in the direction k2
    inline static poly prediff(const KEY& k1, const LET& k2)
    {
        typename KEY::iterator it;
        KEY k(k1);
        it = k.find(k2);
        poly result;// zero
        if (it != k.end()) {
            if (it->second == 1) {
                k.erase(it);
                result = poly(k);
            }
            else {
                SCA coeff = SCA((it->second)--);
                poly temp1(k);
                poly temp2(coeff);
                result = temp2 * temp1;
            }
        }
        return result;
    }

public:
    /// Partial differentiation of a polynomail in the direction k2.
    inline static poly diff(const poly& p1, const LET& k2)
    {
        poly result;
        const_iterator it;
        for (it = p1.begin(); it != p1.end(); ++it) {
            result += poly(it->value()) * prediff(it->key(), k2);
        }
        return result;
    }

    /// Computes the truncated exponential of arg

    /// The result is 1 + arg + arg^2/2! + ... + arg^n/n! where n = max_degree
    inline friend poly exp(const poly& arg, DEG max_degree = 3)
    {
        static KEY kunit;
        poly result(kunit);
        for (DEG i = max_degree; i >= 1; --i) {
            result.mul_scal_div(arg, (RAT)i);
            result += (poly)kunit;
        }
        return result;
    }

    /// Computes the truncated logarithm of a poly instance.

    /// Computes the truncated log of arg up to degree max_degree
    /// The coef. of the constant term (empty word in the monoid) of arg
    /// is forced to 1.
    /// log(arg) = log(1+x) = x - x^2/2 + ... + (-1)^(n+1) x^n/n.
    /// max_degree must be > 0.
    inline friend poly log(const poly& arg, DEG max_degree = 3)
    {
        static KEY kunit;
        poly tunit(kunit);
        poly x(arg);
        iterator it = x.find(kunit);
        if (it != x.end()) {
            x.erase(it);
        }
        poly result;
        for (DEG i = max_degree; i >= 1; --i) {
            if (i % 2 == 0) {
                result.sub_scal_div(tunit, (RAT)i);
            }
            else {
                result.add_scal_div(tunit, (RAT)i);
            }
            result *= x;
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
// Include once wrapper
#endif// DJC_COROPA_LIBALGEBRA_POLYNOMIALSH_SEEN

// EOF.
