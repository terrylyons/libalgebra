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


struct one_method_multiplication_tag
{
};
struct two_method_multiplication_tag
{
};

template<typename SelectBasis>
struct basis_multiplication_selector
{
    typedef one_method_multiplication_tag tag;

    template<typename Basis, typename Coeffs, typename Transform>
    struct key_operator
    {

        typedef typename Basis::KEY KEY;
        typedef typename Coeffs::S S;

        /// Trivial constructor
        key_operator() : m_transform()
        {}

        /// Passthrough constructor for transform
        template<typename Arg>
        key_operator(Arg a) : m_transform(a)
        {}

        template<typename Vector>
        inline void operator()(
                Vector &result,
                const KEY &lhs_key,
                const S &lhs_val,
                const KEY &rhs_key,
                const S &rhs_val
        )
        {
            result.add_scal_prod(
                    Vector::basis.prod(lhs_key, rhs_key),
                    m_transform(lhs_val * rhs_val)
            );
        }

    private:
        Transform m_transform;
    };
};





/// A class to store and manipulate associative algebras elements.
/**
The template class BASIS must
(1) Satisfies the assumptions made by the sparse_vector template class.
(2) Provides two member functions
DEG BASIS::degree(const KEY&) const
BASIS::prod(const KEY&, const KEY&) with a return type suitable for
use as the first arg of sparse_vector::add_scal_prod(); it can be a key or a sparse vector for example
(3) The sparse_vector::MAP class must provide the swap() member function.
*/
template<typename Basis, typename Coeff, typename VectorType>
class algebra : public vectors::vector<Basis, Coeff, VectorType>
{
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
    using VECT::coefficient_field;


    static const DEG MAX_DEGREE = BASIS::MAX_DEGREE;




    // Transformations

    /// function objects for doing nothing to a scalar
    struct scalar_passthrough
    {
        SCALAR operator()(const SCALAR &arg)
        {
            return arg;
        }
    };

    /// function object for changing the sign of a scalar
    struct scalar_minus
    {
        SCALAR operator()(const SCALAR &arg)
        {
            return -arg;
        }
    };

    /// function object for post-multiplying a scalar by a scalar
    struct scalar_post_mult
    {
    private:
        SCALAR mFactor;
    public:
        scalar_post_mult(const SCALAR Factor = VECT::one)
                : mFactor(Factor)
        {}

        SCALAR operator()(const SCALAR arg)
        {
            return arg * mFactor;
        }
    };

    /// function object for post-multiplying a scalar by stored version of the scalar 1 / rational
    struct rational_post_div
    {
    private:
        SCALAR mFactor;
    public:
        rational_post_div(const RATIONAL Factor = VECT::one)
                : mFactor(VECT::one / Factor)
        {}

        SCALAR operator()(const SCALAR arg)
        {
            return arg * mFactor;
        }
    };

    // used to index types
    template<unsigned int T>
    struct identity
    {
    };


    template<typename Transform>
    inline void buffered_apply_binary_transform(
            algebra &result,
            const algebra &rhs,
            Transform fn,
            one_method_multiplication_tag
    ) const
    {
        typename basis_multiplication_selector<BASIS>::
        template key_operator<BASIS, Coeff, Transform>
                key_fn(fn);
        VECT::buffered_apply_binary_transform(result, rhs, key_fn);
    }

    template<typename Transform>
    inline void buffered_apply_binary_transform(
            algebra &result,
            const algebra &rhs,
            Transform fn,
            two_method_multiplication_tag
    ) const
    {
        typename basis_multiplication_selector<BASIS>::
        template index_operator<BASIS, Coeff, Transform> index_fn(fn);
        typename basis_multiplication_selector<BASIS>::
        template key_operator<BASIS, Coeff, Transform> key_fn(fn);
        VECT::buffered_apply_binary_transform(result, rhs, key_fn, index_fn);
    }


    /// multiplies *this and rhs adding it to result
    template<unsigned DEPTH1>
    inline void bufferedmultiplyandadd(const algebra &rhs, algebra &result) const
    {
        typename basis_multiplication_selector<BASIS>::tag tag;
        buffered_apply_binary_transform(result, rhs, scalar_passthrough(), tag);
    }

public:
    /// multiplies *this and rhs subtracting it from result
    template<unsigned DEPTH1>
    inline void bufferedmultiplyandsub(const algebra &rhs, algebra &result) const
    {
        typename basis_multiplication_selector<BASIS>::tag tag;
        buffered_apply_binary_transform(result, rhs, scalar_minus(), tag);
    }

public:
    struct wrapscalar
    {
        const SCALAR &hidden;

        wrapscalar(const SCALAR &s)
                : hidden(s)
        {}
    };

    struct wraprational
    {
        const RATIONAL &hidden;

        wraprational(const RATIONAL &s)
                : hidden(s)
        {}
    };

    /// multiplies  *this and rhs adds it * s to result
    template<unsigned DEPTH1>
    inline void bufferedmultiplyandsmult(const algebra &rhs, const wrapscalar &ss, algebra &result) const
    {
        typename basis_multiplication_selector<BASIS>::tag tag;
        buffered_apply_binary_transform(result, rhs, scalar_post_mult(ss.hidden), tag);
    }

public:
    /// multiplies  *this and rhs adds it * s to result
    template<unsigned DEPTH1>
    inline void bufferedmultiplyandsdiv(const algebra &rhs, const wraprational &ss, algebra &result) const
    {
        typename basis_multiplication_selector<BASIS>::tag tag;
        buffered_apply_binary_transform(result, rhs, rational_post_div(ss.hidden), tag);
    }

public:
    /// Default constructor.
    /**
    Constructs an empty algebra element.
    */
    algebra(void)
    {}

    /// Copy constructor.
    algebra(const algebra &a)
            : VECT(a)
    {}

    /// Constructs an algebra instance from a sparse_vector.
    algebra(const VECT &v)
            : VECT(v)
    {}

    /// Unidimensional constructor.
    explicit algebra(const KEY &k, const SCALAR &s = VECT::one)
            : VECT(k, s)
    {}

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
        algebra result;
            bufferedmultiplyandadd<MAX_DEGREE>(rhs, result);
        this->swap(result);
        return *this;
    }

    /// Binary version of the product of algebra instances.
    inline __DECLARE_BINARY_OPERATOR(algebra, *, *=, algebra);

    /// Adds to the instance a product of algebra instances.
    inline algebra &add_mul(const algebra &a, const algebra &b)
    {
        a.bufferedmultiplyandadd<MAX_DEGREE>(b, *this);
        return *this;
    }

    /// Subtracts to the instance a product of algebra instances.
    inline algebra &sub_mul(const algebra &a, const algebra &b)
    {
        a.bufferedmultiplyandsub<MAX_DEGREE>(b, *this);
        return *this;
    }

    /// Multiplies the instance by (algebra instance)*s.
    inline algebra &mul_scal_prod(const algebra &rhs, const SCALAR &s)
    {
        algebra result;
            bufferedmultiplyandsmult<MAX_DEGREE>(rhs, wrapscalar(s), result);
        this->swap(result);
        return *this;
    }

    algebra& mul_scal_prod(const algebra& rhs, const RATIONAL&s, const DEG depth)
    {
        algebra result;
        buffered_apply_binary_transform(result, rhs, scalar_post_mult(s), depth);
        this->swap(result);
        return *this;
    }

    /// Multiplies the instance by (algebra instance)/s.
    inline algebra &mul_scal_div(const algebra &rhs, const RATIONAL &s)
    {
        algebra result;
            bufferedmultiplyandsdiv<MAX_DEGREE>(rhs, wraprational(s), result);
        this->swap(result);
        return *this;
    }

    algebra& mul_scal_div(const algebra& rhs, const RATIONAL&s, const DEG depth)
    {
        algebra result;
        buffered_apply_binary_transform(result, rhs, rational_post_div(s), depth);
        this->swap(result);
        return *this;
    }


    /// Returns an instance of the commutator of two algebra instances.
    inline friend algebra commutator(const algebra &a, const algebra &b)
    { // Returns a * b - b * a
        algebra result;
        a.bufferedmultiplyandadd<MAX_DEGREE>(b, result);
        b.bufferedmultiplyandsub<MAX_DEGREE>(a, result);
        return result;
    }

    /// Returns a truncated version of the instance, by using basis::degree().
    inline algebra truncate(const DEG min, const DEG max) const
    {
        algebra result;
        const_iterator i;
        for (i = VECT::begin(); i != VECT::end(); ++i)
            if ((VECT::basis.degree(i->key()) >= min) && (VECT::basis.degree(i->key()) <= max))
                result[i->key()] = i->value();
        return result;
    }

    /// Returns the degree of the instance by using basis:degree()
    inline DEG degree(void) const
    {
        DEG result(0);
        const_iterator i;
        for (i = VECT::begin(); i != VECT::end(); ++i)
            result = std::max(result, VECT::basis.degree(i->key()));
        return result;
    }
};

// Include once wrapper
#endif // DJC_COROPA_LIBALGEBRA_ALGEBRAH_SEEN

//EOF.
