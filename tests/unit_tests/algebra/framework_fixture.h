//
// Created by sam on 12/02/2021.
//

#ifdef ALGEBRA_TESTS_VECT_TYPE

struct Fixture
{
    typedef alg_types<5, 5, Rational> ALG_TYPES;

    typedef typename ALG_TYPES::COEFF field;

    typedef SimpleIntegerBasis<5, field> BASIS;
    typedef alg::vectors::vector<BASIS, field,
    ALGEBRA_TESTS_VECT_TYPE> VECT;
    typedef alg::algebra<BASIS, field,
    pointwise_multiplication<BASIS>,
    ALGEBRA_TESTS_VECT_TYPE> ALG;

    typedef typename BASIS::KEY KEY;
    typedef typename field::S S;
    typedef typename field::Q Q;

    typedef mt19937 RNG;

    RNG m_rng;

    S rand_scalar(S lower_bound = S(-5), S upper_bound = S(5))
    {
        RNG::result_type numerator = m_rng(), denominator = m_rng.max();

        S uniform01 = S(numerator) / S(denominator);
        return (upper_bound - lower_bound)*uniform01 + lower_bound;
    }

    ALG rand_vec(S lower_bound = S(-5), S upper_bound = S(5))
    {
        ALG v;
        for (KEY i=0; i<BASIS::dimension; ++i)
            v[i] = rand_scalar(lower_bound, upper_bound);
        return v;
    }

    KEY rand_key() {
        return (m_rng() % BASIS::dimension);
    }


};
#endif
