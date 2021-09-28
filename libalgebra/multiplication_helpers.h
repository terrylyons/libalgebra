//
// Created by sam on 05/07/2021.
//

#ifndef LIBALGEBRA_MULTIPLICATION_HELPERS_H
#define LIBALGEBRA_MULTIPLICATION_HELPERS_H

namespace alg {
namespace mult {

struct scalar_passthrough
{
    template <typename S> LA_CONSTEXPR S operator()(S val) const { return val; }
};

template <typename Coeff> struct scalar_minus
{
    typedef typename Coeff::SCA scalar_t;

    LA_CONSTEXPR scalar_t operator()(scalar_t val) const
    {
        return Coeff::uminus(val);
    }
};

template <typename Coeff> struct scalar_post_mult
{
    typedef typename Coeff::SCA scalar_t;

    LA_CONSTEXPR explicit scalar_post_mult(const scalar_t factor = Coeff::one) : m_factor(factor) {}

    LA_CONSTEXPR scalar_t operator()(scalar_t val) const
    {
        return Coeff::mul(val, m_factor);
    }

private:
    scalar_t m_factor;
};

template <typename Coeff> struct rational_post_div
{
    typedef typename Coeff::SCA scalar_t;
    typedef typename Coeff::RAT rational_t;

    LA_CONSTEXPR explicit rational_post_div(const rational_t factor) : m_factor(Coeff::div(Coeff::one, factor)) {}

    LA_CONSTEXPR scalar_t operator()(scalar_t val) const
    {
        return Coeff::mul(val, m_factor);
    }

private:
    scalar_t m_factor;
};

} // namespace mult
} // namespace alg

#endif // LIBALGEBRA_MULTIPLICATION_HELPERS_H
