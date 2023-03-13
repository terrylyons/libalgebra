//
// Created by sam on 03/02/2021.
//
#if defined(_VECTOR_TYPE)

struct Fixture {
    typedef alg_types<2, 2, Rational> FIELD;
    typedef SimpleIntegerBasis<5, typename FIELD::RAT> BASIS;
    typedef _VECTOR_TYPE<BASIS
#ifdef LIBALGEBRA_VECTORS_H
                         ,
                         typename FIELD::COEFF
#endif
#ifdef _VECTOR_TYPE_ADDITIONAL_PARAMS
                         ,
                         _VECTOR_TYPE_ADDITIONAL_PARAMS
#endif
                         >
            VECT;

    typedef typename FIELD::S S;
    typedef typename FIELD::Q Q;
    typedef typename BASIS::KEY KEY;
    typedef typename FIELD::LET LET;

    using scalar_t = S;
    using IDIMN = std::ptrdiff_t;

    typedef alg::free_tensor_basis<2, 5> TBASIS;
    typedef typename TBASIS::KEY TKEY;
    typedef _VECTOR_TYPE<TBASIS
#ifdef LIBALGEBRA_VECTORS_H
                         ,
                         typename FIELD::COEFF
#endif
#ifdef _TVECTOR_TYPE_ADDITIONAL_PARAMS
                         ,
                         _TVECTOR_TYPE_ADDITIONAL_PARAMS
#endif
                         >
            TVECT;

    typedef mt19937 RNG;
    RNG m_rng;

    template<typename Transform>
    class ft_index_operator
    {
        Transform m_transform;

    public:
        ft_index_operator(Transform t)
            : m_transform(t) {}

        void operator()(scalar_t* result_ptr, scalar_t const* lhs_ptr, scalar_t const* rhs_ptr, DIMN const lhs_target,
                        DIMN const rhs_target, bool assign = false)
        {
            scalar_t lhs;
            if (assign) {
                for (DIMN i = 0; i < lhs_target; ++i) {
                    lhs = lhs_ptr[i];
                    for (DIMN j = 0; j < rhs_target; ++j) {
                        *(result_ptr++) = m_transform(lhs * rhs_ptr[j]);
                    }
                }
            }
            else {
                for (DIMN i = 0; i < lhs_target; ++i) {
                    lhs = lhs_ptr[i];
                    for (DIMN j = 0; j < rhs_target; ++j) {
                        *(result_ptr++) += m_transform(lhs * rhs_ptr[j]);
                    }
                }
            }
        }
    };

    template<typename Transform>
    class ft_key_operator
    {
        Transform m_transform;

    public:
        ft_key_operator(Transform t)
            : m_transform(t) {}

        template<typename Vector>
        void
        operator()(Vector& result, typename Vector::KEY const& lhs_key, scalar_t const& lhs_val,
                   typename Vector::KEY const& rhs_key, scalar_t const& rhs_val)
        {
            result.add_scal_prod(lhs_key * rhs_key, m_transform(lhs_val * rhs_val));
        }
    };

    S rand_scalar(S lower_bound = S(-5), S upper_bound = S(5))
    {
        RNG::result_type numerator = m_rng(), denominator = m_rng.max();

        S uniform01 = S(numerator) / S(denominator);
        return (upper_bound - lower_bound) * uniform01 + lower_bound;
    }

    VECT rand_vec(S lower_bound = S(-5), S upper_bound = S(5))
    {
        VECT v;
        for (KEY i = 0; i < BASIS::dimension; ++i)
            v[i] = rand_scalar(lower_bound, upper_bound);
        return v;
    }

    TVECT rand_tvec(const DEG deg, S lower_bound = S(-5), S upper_bound = S(5))
    {
        TVECT v;
        TKEY k;
        TBASIS& basis(v.basis);

        while (basis.degree(k) <= deg) {
            v[k] = rand_scalar(lower_bound, upper_bound);
            k = basis.nextkey(k);
        }

        return v;
    }

    KEY rand_key()
    {
        return (m_rng() % BASIS::dimension);
    }
};

#endif