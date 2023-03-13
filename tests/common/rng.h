//
// Created by sam on 24/02/2021.
//

#ifndef LIBALGEBRAUNITTESTS_RNG_H
#define LIBALGEBRAUNITTESTS_RNG_H

#include <boost/random.hpp>
#include <limits>
#include <random>

#include <libalgebra/rational_coefficients.h>

#if __cplusplus >= 201103L

typedef std::mt19937 mt19937;
#define NORMAL_DIST std::normal_distribution
#define UNIFORM_INT_DIST std::uniform_int_distribution
#else

typedef boost::mt19937 mt19937;
#define NORMAL_DIST boost::normal_distribution
#define UNIFORM_INT_DIST boost::random::uniform_int_distribution
#endif

namespace la_testing {

constexpr unsigned rebase_const = std::numeric_limits<unsigned>::max();

template<typename Rational>
struct uniform_rational_distribution {

    using integer_t = decltype(boost::multiprecision::numerator(std::declval<Rational>()));
    using integer_dist_t = boost::random::uniform_int_distribution<integer_t>;

    uniform_rational_distribution(Rational min, Rational max)
        : m_denom(static_cast<integer_t>(rebase_const) * static_cast<integer_t>(rebase_const) * boost::multiprecision::denominator(min) * boost::multiprecision::denominator(max)),
          m_int_dist(boost::multiprecision::numerator(min) * (static_cast<integer_t>(rebase_const) * boost::multiprecision::denominator(max)),
                     boost::multiprecision::numerator(max) * (static_cast<integer_t>(rebase_const) * boost::multiprecision::denominator(min)))
    {}

    template<typename Rng>
    Rational operator()(Rng& rng) const
    {
        integer_t numerator = m_int_dist(rng);
        return Rational(numerator) / m_denom;
    }

private:
    integer_t m_denom;
    integer_dist_t m_int_dist;
};

}// namespace la_testing

#endif//LIBALGEBRAUNITTESTS_RNG_H
