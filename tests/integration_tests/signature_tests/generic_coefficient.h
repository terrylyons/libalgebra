//
// Created by sam on 11/03/2021.
//

#ifndef LIBALGEBRAUNITTESTS_GENERIC_COEFFICIENT_H
#define LIBALGEBRAUNITTESTS_GENERIC_COEFFICIENT_H

#if __cplusplus >= 201103UL
#define CONSTEXPR constexpr
#else
#define CONSTEXPR
#endif

#include "tests/common/rng.h"

template<typename Integer>
struct generic_coefficient {
    Integer numerator;
    Integer denominator;

    typedef typename alg_types<2, 2, Rational>::SCA rational_t;

    generic_coefficient() : numerator(0), denominator(1)
    {}

    generic_coefficient(Integer a, Integer b) : numerator(a), denominator(b)
    {}

    explicit generic_coefficient(Integer i) : numerator(i), denominator(1)
    {}

    generic_coefficient& operator=(const generic_coefficient& other)
    {
        numerator = other.numerator;
        denominator = other.denominator;
        return *this;
    }

    explicit operator float() const
    {
        return float(numerator) / float(denominator);
    }

    explicit operator double() const
    {
        return double(numerator) / float(denominator);
    }

    explicit operator rational_t() const
    {
        return rational_t(numerator) / rational_t(denominator);
    }
};

#endif//LIBALGEBRAUNITTESTS_GENERIC_COEFFICIENT_H
