//
// Created by sam on 05/07/2021.
//

#ifndef LIBALGEBRA_MULTIPLICATION_HELPERS_H
#define LIBALGEBRA_MULTIPLICATION_HELPERS_H

namespace alg {
namespace mult {

/// Identity function on scalars
struct scalar_passthrough {
    template<typename S>
    constexpr S operator()(S val) const noexcept
    {
        return val;
    }
};

/// Scalar minus function on scalars
template<typename Coeff>
struct scalar_minus {
    typedef typename Coeff::SCA scalar_t;

    constexpr scalar_t operator()(scalar_t val) const noexcept(noexcept(Coeff::uminus(val)))
    {
        return Coeff::uminus(val);
    }
};

/// Post multiply scalar by fixed scalar
template<typename Coeff>
struct scalar_post_mult {
    typedef typename Coeff::SCA scalar_t;

    constexpr explicit scalar_post_mult(const scalar_t factor = Coeff::one) noexcept(noexcept(scalar_t(factor)))
        : m_factor(factor)
    {}

    constexpr scalar_t operator()(scalar_t val) const noexcept(noexcept(Coeff::mul(val, std::declval<scalar_t>())))
    {
        return Coeff::mul(val, m_factor);
    }

private:
    scalar_t m_factor;
};

/// Post divide scalar by fixed rational
template<typename Coeff>
struct rational_post_div {
    typedef typename Coeff::SCA scalar_t;
    typedef typename Coeff::RAT rational_t;

    constexpr explicit rational_post_div(const rational_t factor)
        : m_factor(Coeff::div(Coeff::one, factor))
    {}

    constexpr scalar_t operator()(scalar_t val) const noexcept(noexcept(Coeff::mul(val, std::declval<scalar_t>())))
    {
        return Coeff::mul(val, m_factor);
    }

private:
    scalar_t m_factor;
};

}// namespace mult
}// namespace alg

#endif// LIBALGEBRA_MULTIPLICATION_HELPERS_H
