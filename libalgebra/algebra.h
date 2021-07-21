/* *************************************************************

Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurk� and Arend Janssen.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */

//  algebra.h

// Include once wrapper
#ifndef DJC_COROPA_LIBALGEBRA_ALGEBRAH_SEEN
#define DJC_COROPA_LIBALGEBRA_ALGEBRAH_SEEN

/// Temporary implementation of a basis-level key_transform for multiplication.

struct one_method_multiplication_tag {};
struct two_method_multiplication_tag {};

template <typename SelectBasis> struct basis_multiplication_selector
{
    typedef one_method_multiplication_tag tag;

    template <typename Basis, typename Coeffs, typename Transform> struct key_operator
    {

        typedef typename Basis::KEY KEY;
        typedef typename Coeffs::S S;

        /// Trivial constructor
        key_operator() : m_transform() {}

        /// Passthrough constructor for transform
        template <typename Arg> key_operator(Arg a) : m_transform(a) {}

        template <typename Vector> inline void
        operator()(Vector &result, const KEY &lhs_key, const S &lhs_val, const KEY &rhs_key, const S &rhs_val)
        {
            result.add_scal_prod(Vector::basis.prod(lhs_key, rhs_key), m_transform(lhs_val * rhs_val));
        }

    private:
        Transform m_transform;
    };
};

/* The algebra is a vector with an additional multiplication operation which is
 * "compatible" with the addition and scalar multiplication. We attempt to mimic
 * this with the algebra template class below. For this, we need to define
 * "multiplication" types for each of the algebras that we wish to implement.
 * Here is a simple example of what a multiplication class should look like.
 *
 * template <typename Coeff>
 * class multiplication {
 *
 *      typedef typename Coeff::SCA scalar_t;
 *
 *      template <typename Transform>
 *      class index_operator {
 *          Transform m_transform;
 *
 *      public:
 *
 *          index_operator(Transform t) : m_transform(t) {}
 *
 *
 *          void operator()(
 *              scalar_t *result_ptr,
 *              scalar_t const* lhs_ptr,
 *              scalar_t const* rhs_ptr,
 *              DIMN const lhs_target,
 *              DIMN const rhs_target,
 *              bool assign = false
 *          );
 *      };
 *
 *      template <typename Transform>
 *      class key_operator {
 *          Transform m_transform;
 *
 *          typedef ...    key_t;
 *
 *      public:
 *          key_transform(Transform t) : m_transform(t) {}
 *
 *          template <typename Vector>
 *          void operator()(
 *              Vector& result,
 *              key_t const& lhs_key,
 *              scalar_t const& lhs_val,
 *              key_t const& rhs_key,
 *              scalar_t const& rhs_val
 *          );
 *      };
 *
 *
 * public:
 *
 *      template <typename Algebra, typename Operator>
 *      Algebra& multiply_and_add(Algebra& result, Algebra const& lhs, Algebra
 * const& rhs, Operator op) const { key_transform<scalar_t, Operator> kt(op);
 *          index_transform<scalar_t, Operator> it(op);
 *          lhs.buffered_apply_binary_transform(result, rhs, kt, it);
 *          return result;
 *      }
 *
 *      template <typename Algebra, typename Operator>
 *      Algebra multiply(Algebra const& lhs, Algebra const& rhs, Operator op)
 * const { typedef typename Algebra::SCALAR scalar_t; Algebra result;
 *          multiply_and_add(result, lhs, rhs, op);
 *          return result;
 *      }
 *
 *      template <typename Algebra, typename Operator>
 *      Algebra& multiply_inplace(Algebra& lhs, Algebra const& rhs, Operator op)
 * const { typedef typename Algebra::SCALAR scalar_t; key_transform<scalar_t,
 * Operator> kt(op); index_transform<scalar_t, Operator> it(op);
 *          lhs.unbuffered_apply_binary_transform(rhs, kt, it);
 *          return lhs;
 *      }
 * };
 *
 *
 *
 *
 *
 */

/// A class to store and manipulate associative algebras elements.
/**
The template class BASIS must
(1) Satisfies the assumptions made by the sparse_vector template class.
(2) Provides two member functions
DEG BASIS::degree(const KEY&) const
BASIS::prod(const KEY&, const KEY&) with a return type suitable for
use as the first arg of sparse_vector::add_scal_prod(); it can be a key or a
sparse vector for example (3) The sparse_vector::MAP class must provide the
swap() member function.
*/
template <typename Basis, typename Coeff, typename Multiplication, typename VectorType>
class algebra : public vectors::vector<Basis, Coeff, VectorType>
{

    typedef mult::scalar_passthrough scalar_passthrough;
    typedef mult::scalar_minus <Coeff> scalar_minus;
    typedef mult::scalar_post_mult <Coeff> scalar_post_mult;
    typedef mult::rational_post_div <Coeff> rational_post_div;

public:
    typedef Basis BASIS;
    /// The inherited sparse vector type.
    typedef vectors::vector <Basis, Coeff, VectorType> VECT;
    /// Import of the iterator type from sparse_vector.
    typedef typename VECT::iterator iterator;
    /// Import of the constant iterator type from sparse_vector.
    typedef typename VECT::const_iterator const_iterator;
    /// Import of the KEY type from sparse_vector.
    typedef typename VECT::KEY KEY;
    /// Import of the SCALAR type from sparse_vector.
    typedef typename VECT::SCALAR SCALAR;
    /// Import of the RATIONAL type from sparse_vector.
    typedef typename VECT::RATIONAL RATIONAL;
    using typename VECT::coefficient_field;

    typedef Multiplication multiplication_t;

    static const DEG MAX_DEGREE = BASIS::MAX_DEGREE;

private:
    static multiplication_t s_multiplication;


public:
    /// Default constructor.
    /**
    Constructs an empty algebra element.
    */
    algebra(void) {}

    /// Copy constructor.
    algebra(const algebra &a) : VECT(a) {}

    /// Constructs an algebra instance from a sparse_vector.
    algebra(const VECT &v) : VECT(v) {}

    /// Unidimensional constructor.
    explicit algebra(const KEY &k, const SCALAR &s = VECT::one) : VECT(k, s) {}

public:
    /// Multiplies the instance with scalar s.
    inline algebra &operator*=(const SCALAR &s)
    {
        VECT::operator*=(s);
        return *this;
    }

    /// Divides the instance by scalar s.
    inline algebra &operator/=(const RATIONAL &s)
    {
        VECT::operator/=(s);
        return *this;
    }

    /// Ensures that the return type is an instance of algebra.
    inline __DECLARE_BINARY_OPERATOR(algebra, *, *=, SCALAR)

    /// Ensures that the return type is an instance of algebra.
    inline __DECLARE_BINARY_OPERATOR(algebra, /, /=, SCALAR)

    /// Ensures that the return type is an instance of algebra.
    inline __DECLARE_BINARY_OPERATOR(algebra, +, +=, algebra)

    /// Ensures that the return type is an instance of algebra.
    inline __DECLARE_BINARY_OPERATOR(algebra, -, -=, algebra)

    /// Ensures that the return type is an instance of algebra.
    inline __DECLARE_UNARY_OPERATOR(algebra, -, -, VECT);

    /// Multiplies the instance by an instance of algebra.
    inline algebra &operator*=(const algebra &rhs)
    {
        return s_multiplication.multiply_inplace(*this, rhs, scalar_passthrough());
    }

    /// Binary version of the product of algebra instances.
    // inline __DECLARE_BINARY_OPERATOR(algebra, *, *=, algebra);
    algebra operator*(algebra const &rhs) const
    {
        return s_multiplication.multiply(*this, rhs, scalar_passthrough());
    }

    /// Adds to the instance a product of algebra instances.
    inline algebra &add_mul(const algebra &a, const algebra &b)
    {
        s_multiplication.multiply_and_add(*this, a, b, scalar_passthrough());
        return *this;
    }

    /// Subtracts to the instance a product of algebra instances.
    inline algebra &sub_mul(const algebra &a, const algebra &b)
    {
        return s_multiplication.multiply_and_add(*this, a, b, scalar_minus());
    }

    /// Multiplies the instance by (algebra instance)*s.
    inline algebra &mul_scal_prod(const algebra &rhs, const SCALAR &s)
    {
        return s_multiplication.multiply_inplace(*this, rhs, scalar_post_mult(s));
    }

    algebra &mul_scal_prod(const algebra &rhs, const RATIONAL &s, const DEG depth)
    {
        return s_multiplication.multiply_inplace(*this, rhs, rational_post_div(s), depth);
    }

    /// Multiplies the instance by (algebra instance)/s.
    inline algebra &mul_scal_div(const algebra &rhs, const RATIONAL &s)
    {
        return s_multiplication.multiply_inplace(*this, rhs, rational_post_div(s));
    }

    algebra &mul_scal_div(const algebra &rhs, const RATIONAL &s, const DEG depth)
    {
        return s_multiplication.multiply_inplace(*this, rhs, rational_post_div(s), depth);
    }

    /// Returns an instance of the commutator of two algebra instances.
    inline friend algebra commutator(const algebra &a, const algebra &b)
    { // Returns a * b - b * a
        algebra result;
        s_multiplication.multiply_and_add(result, a, b, scalar_passthrough());
        return s_multiplication.multiply_and_add(result, b, a, scalar_minus());
    }

    /// Returns a truncated version of the instance, by using basis::degree().
    inline algebra truncate(const DEG min, const DEG max) const
    {
        algebra result;
        const_iterator i;
        for (i = VECT::begin(); i != VECT::end(); ++i) {
            if ((VECT::basis.degree(i->key()) >= min) && (VECT::basis.degree(i->key()) <= max)) {
                result[i->key()] = i->value();
            }
        }
        return result;
    }


};

template <typename B, typename C, typename M, typename V> M algebra<B, C, M, V>::s_multiplication;

// Include once wrapper
#endif // DJC_COROPA_LIBALGEBRA_ALGEBRAH_SEEN

// EOF.
