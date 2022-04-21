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

namespace alg {

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

/**
 * @brief Class to store and manipulate associative algebra elements
 *
 * An algebra is a vector space that is given a compatible multiplication operation.
 * An algebra class is constructed in the same way. We provide the basis and coefficients
 * required to instantiate a vector class and provide an additional multiplication class
 * that implements the multiplication via key operators (acting on sparse data) and index
 * operators (operating on dense data).
 *
 * @tparam Basis Basis of underlying vector space
 * @tparam Coeff Coefficient field of underlying vector space
 * @tparam Multiplication Multiplication operation
 * @tparam VectorType Underlying vector data type to use; e.g. dense_vector or sparse_vector
 */
template<typename Basis, typename Coeff, typename Multiplication,
         template<typename, typename, typename...> class VectorType = alg::vectors::template_vector_type_selector<Basis, Coeff>::template type,
         typename... Args>
class algebra : public vectors::vector<Basis, Coeff, VectorType, Args...>
{

    typedef mult::scalar_passthrough scalar_passthrough;
    typedef mult::scalar_minus<Coeff> scalar_minus;
    typedef mult::scalar_post_mult<Coeff> scalar_post_mult;
    typedef mult::rational_post_div<Coeff> rational_post_div;

public:
    typedef Basis BASIS;
    /// The inherited sparse vector type.
    typedef vectors::vector<Basis, Coeff, VectorType, Args...> VECT;
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

private:
    static multiplication_t s_multiplication;

public:
    /// Default constructor.
    /**
    Constructs an empty algebra element.
    */
    algebra()
        : VECT()
    {}

    /// Copy constructor.
    algebra(const algebra& a)
        : VECT(a)
    {}

    /// Move constructor
    algebra(algebra&& a) noexcept
        : VECT(std::move(a))
    {}

    /// Constructs an algebra instance from a sparse_vector.
    explicit algebra(const VECT& v)
        : VECT(v)
    {}

    /// Unidimensional constructor.
    explicit algebra(const KEY& k, const SCALAR& s = VECT::one)
        : VECT(k, s)
    {}

    /// Constructor from pointers to data range
    algebra(SCALAR const* begin, SCALAR const* end)
        : VECT(begin, end)
    {}

    /// Constructor from pointer to data range with offset
    algebra(DIMN offset, SCALAR const* begin, SCALAR const* end)
        : VECT(offset, begin, end)
    {}

    /// Constructor from pointer to data range with offset
    algebra(DIMN offset, SCALAR* begin, SCALAR* end)
        : VECT(offset, begin, end)
    {}

    algebra& operator=(const algebra&) = default;
    algebra& operator=(algebra&&) noexcept = default;

public:
    /// Multiplies the instance with scalar s.
    inline algebra& operator*=(const SCALAR& s)
    {
        VECT::operator*=(s);
        return *this;
    }

    /// Divides the instance by scalar s.
    inline algebra& operator/=(const RATIONAL& s)
    {
        VECT::operator/=(s);
        return *this;
    }

    /// Ensures that the return type is an instance of algebra.
    inline algebra operator*(const SCALAR& rhs) const
    {
        algebra result(*this);
        result *= rhs;
        return result;
    }

    /// Ensures that the return type is an instance of algebra.
    inline algebra operator/(const SCALAR& rhs) const
    {
        algebra result(*this);
        result /= rhs;
        return result;
    }

    /// Ensures that the return type is an instance of algebra.
    inline algebra operator+(const algebra& rhs) const
    {
        algebra result(*this);
        result += rhs;
        return result;
    }

    /// Ensures that the return type is an instance of algebra.
    inline algebra operator-(const algebra& rhs) const
    {
        algebra result(*this);
        result -= rhs;
        return result;
    }

    /// Ensures that the return type is an instance of algebra.
    inline algebra operator-() const { return algebra(VECT::operator-()); };

    /// Multiplies the instance by an instance of algebra.
    inline algebra& operator*=(const algebra& rhs)
    {
        return s_multiplication.multiply_inplace(*this, rhs, scalar_passthrough());
    }

    /// Binary version of the product of algebra instances.
    // inline __DECLARE_BINARY_OPERATOR(algebra, *, *=, algebra);
    algebra operator*(algebra const& rhs) const
    {
        return s_multiplication.multiply(*this, rhs, scalar_passthrough());
    }

    /// Adds to the instance a product of algebra instances.
    inline algebra& add_mul(const algebra& a, const algebra& b)
    {
        s_multiplication.multiply_and_add(*this, a, b, scalar_passthrough());
        return *this;
    }

    /// Subtracts to the instance a product of algebra instances.
    inline algebra& sub_mul(const algebra& a, const algebra& b)
    {
        return s_multiplication.multiply_and_add(*this, a, b, scalar_minus());
    }

    /// Multiplies the instance by (algebra instance)*s.
    inline algebra& mul_scal_prod(const algebra& rhs, const SCALAR& s)
    {
        return s_multiplication.multiply_inplace(*this, rhs, scalar_post_mult(s));
    }

    /// Multiplies the instance by (algebra instance)*s up to maximum depth
    algebra& mul_scal_prod(const algebra& rhs, const RATIONAL& s, const DEG depth)
    {
        return s_multiplication.multiply_inplace(*this, rhs, rational_post_div(s), depth);
    }

    /// Multiplies the instance by (algebra instance)/s.
    inline algebra& mul_scal_div(const algebra& rhs, const RATIONAL& s)
    {
        return s_multiplication.multiply_inplace(*this, rhs, rational_post_div(s));
    }

    /// Multiplies the instance by (algebra instance)/s up to maximum depth
    algebra& mul_scal_div(const algebra& rhs, const RATIONAL& s, const DEG depth)
    {
        return s_multiplication.multiply_inplace(*this, rhs, rational_post_div(s), depth);
    }

    /// Returns an instance of the commutator of two algebra instances.
    inline friend algebra commutator(const algebra& a, const algebra& b)
    {// Returns a * b - b * a
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

#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
private:

    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, unsigned int const /* version */) {
        ar & boost::serialization::base_object<VECT>(*this);
    }
#endif
};

template<typename B, typename C, typename M, template<typename, typename, typename...> class V, typename... Args>
M algebra<B, C, M, V, Args...>::s_multiplication;

}// namespace alg
// Include once wrapper
#endif// DJC_COROPA_LIBALGEBRA_ALGEBRAH_SEEN

// EOF.
