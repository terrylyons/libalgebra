/* *************************************************************

Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurkö and Arend Janssen.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */

// multi_polynomial.h

#ifndef multi_polynomialH_SEEN
#define multi_polynomialH_SEEN

namespace alg {

template<typename Coeff>
class multipoly_multiplication
{

    typedef typename Coeff::SCA scalar_t;

    /// The concatenation product of two basis elements.
    /**
    Returns the multi_polynomial obtained by the concatenation product of two keys
    viewed as words of letters. The result is a unidimensional multi_polynomial
    with a unique key (the concatenation of k1 and k2) associated to the +1
    scalar. The already computed products are not stored or remembered.
    */
    template<typename MultiPoly>
    static MultiPoly prod(typename MultiPoly::KEY const& k1, typename MultiPoly::KEY const& k2)
    {

        typedef typename MultiPoly::BASIS basis_t;
        typedef typename MultiPoly::KEY key_t;
        using trait = basis::basis_traits<typename MultiPoly::BASIS>;
        static constexpr DEG max_degree = trait::degree_tag::max_degree;

        MultiPoly result;
        if ((max_degree == 0) || (k1.size() + k2.size() <= max_degree)) {
            key_t concat(k1);
            for (typename key_t::size_type i = 0; i < k2.size(); ++i) {
                concat.push_back(k2[i]);
            }
            result[concat] = Coeff::one;
        }
        return result;
    }


    template<typename MultiPoly, typename Transform>
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
            result.add_scal_prod(prod<MultiPoly>(lhs_key, rhs_key), m_transform(Coeff::mul(lhs_val, rhs_val)));
        }
    };


public:
    template<typename Algebra, typename Operator>
    Algebra& multiply_and_add(Algebra& result, Algebra const& lhs, Algebra const& rhs, Operator op) const
    {
        key_operator<Algebra, Operator> kt(op);
        lhs.buffered_apply_binary_transform(result, rhs, kt);
        return result;
    }

    template<typename Algebra, typename Operator>
    Algebra&
    multiply_and_add(Algebra& result, Algebra const& lhs, Algebra const& rhs, Operator op, DEG const max_depth) const
    {
        key_operator<Algebra, Operator> kt(op);
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
        key_operator<Algebra, Operator> kt(op);
        lhs.unbuffered_apply_binary_transform(rhs, kt);
        return lhs;
    }

    template<typename Algebra, typename Operator>
    Algebra& multiply_inplace(Algebra& lhs, Algebra const& rhs, Operator op, DEG const max_depth) const
    {
        key_operator<Algebra, Operator> kt(op);
        lhs.unbuffered_apply_binary_transform(rhs, kt, max_depth);
        return lhs;
    }
};

/// A specialisation of the algebra class with a free tensor basis.
/**
   Mathematically, the algebra of multi_polynomial instances is a free
   associative algebra. With respect to the inherited algebra class, the
   essential distinguishing feature of this class is the basis class used, and
   in particular the basis::prod() member function. Thus, the most important
   information is in the definition of monomial_basis. Notice that this
   associative algebra of free tensors includes as a sub-algebra the
   associative algebra corresponding to the SCALAR type. This is permitted by
   the existence of empty keys in monomial_basis.
 */
template<typename Coeff, DEG n_letters, DEG max_degree,
         template<typename, typename, typename...> class VectorType,
         typename... Args>
class multi_polynomial : public algebra<
                                 free_monomial_basis<n_letters, max_degree>,
                                 Coeff, multipoly_multiplication<Coeff>,
                                 VectorType, Args...>
{

    typedef multipoly_multiplication<Coeff> multiplication_t;

public:
    /// The basis type.
    typedef free_monomial_basis<n_letters, max_degree> BASIS;
    /// Import of the KEY type.
    typedef typename BASIS::KEY KEY;

    typedef typename Coeff::SCA SCA;
    typedef typename Coeff::RAT RAT;

    /// The algebra type.
    typedef algebra<BASIS, Coeff, multiplication_t, VectorType, Args...> ALG;
    /// The sparse_vector type.
    typedef typename ALG::VECT VECT;
    /// Import of the iterator type.
    typedef typename ALG::iterator iterator;
    /// Import of the constant iterator type.
    typedef typename ALG::const_iterator const_iterator;

public:
    /// Default constructor.
    multi_polynomial()
    {}

    /// Copy constructor.
    multi_polynomial(const multi_polynomial& t)
        : ALG(t)
    {}

    /// Constructs an instance from an algebra instance.
    multi_polynomial(const ALG& a)
        : ALG(a)
    {}

    /// Constructs an instance from a sparse_vector instance.
    multi_polynomial(const VECT& v)
        : ALG(v)
    {}

    /// Constructs a unidimensional instance from a letter and a scalar.
    multi_polynomial(LET letter, const SCA& s) : ALG(VECT::basis.keyofletter(letter), s)
    {
    }

    /// Explicit unidimensional constructor from a given key (basis element).
    explicit multi_polynomial(const KEY& k)
        : ALG(k)
    {}

    /// Explicit unidimensional constructor from a given scalar.
    explicit multi_polynomial(const SCA& s)
        : ALG(VECT::basis.empty_key, s)
    {}

public:
    /// Ensures that the return type is a multi_polynomial.
    inline multi_polynomial operator*(const SCA& rhs) const
    {
        multi_polynomial result(*this);
        result *= rhs;
        return result;
    }

    /// Ensures that the return type is a multi_polynomial.
    inline multi_polynomial
    operator/(const RAT& rhs) const
    {
        multi_polynomial result(*this);
        result /= rhs;
        return result;
    }

    /// Ensures that the return type is a multi_polynomial.
    inline multi_polynomial
    operator*(const multi_polynomial& rhs) const
    {
        multi_polynomial result(*this);
        result *= rhs;
        return result;
    }

    /// Ensures that the return type is a multi_polynomial.
    inline multi_polynomial
    operator+(const multi_polynomial& rhs) const
    {
        multi_polynomial result(*this);
        result += rhs;
        return result;
    }

    /// Ensures that the return type is a multi_polynomial.
    inline multi_polynomial
    operator-(const multi_polynomial& rhs) const
    {
        multi_polynomial result(*this);
        result -= rhs;
        return result;
    }

    /// Ensures that the return type is a multi_polynomial.
    inline multi_polynomial operator-() const { return multi_polynomial(ALG::operator-()); }

    /// Computes the truncated exponential of a multi_polynomial instance.
    inline friend multi_polynomial exp(const multi_polynomial& arg)
    {
        // Computes the truncated exponential of arg
        // 1 + arg + arg^2/2! + ... + arg^n/n! where n = max_degree
        static KEY kunit;
        multi_polynomial result(kunit);
        for (DEG i = max_degree; i >= 1; --i) {
            result.mul_scal_div(arg, (RAT)i);
            result += (multi_polynomial)
                    kunit;
        }
        return result;
    }

    /// Computes the truncated logarithm of a multi_polynomial instance.
    inline friend multi_polynomial log(const multi_polynomial& arg)
    {
        // Computes the truncated log of arg up to degree max_degree
        // The coef. of the constant term (empty word in the monoid) of arg
        // is forced to 1.
        // log(arg) = log(1+x) = x - x^2/2 + ... + (-1)^(n+1) x^n/n.
        // max_degree must be > 0
        static KEY kunit;
        multi_polynomial tunit(kunit);
        multi_polynomial x(arg);
        iterator it = x.find(kunit);
        if (it != x.end()) {
            x.erase(it);
        }
        multi_polynomial result;
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
#endif// DJC_COROPA_LIBALGEBRA_TENSORH_SEEN

// EOF.
