#pragma once

#ifndef LIBALGEBRA_UNITTEST_HELPERS_H
#define LIBALGEBRA_UNITTEST_HELPERS_H

#include <libalgebra/alg_types.h>
#include <boost/type_traits.hpp>


template<bool>
struct check_vec_almost_equal_impl {
    template<class VECT>
    inline static void check(const VECT& expected, const VECT& actual, double tol)
    {
        CHECK_EQUAL(expected, actual);
    }
};

template<>
struct check_vec_almost_equal_impl<false> {
    template<class VECT>
    inline static void check(const VECT& expected, const VECT& actual, double tol)
    {
        VECT err = actual - expected;
        CHECK_CLOSE(0.0, err.NormLInf(), tol);
    }
};

template <class VECT>
inline void check_vec_almost_equal(const VECT& expected, const VECT& actual, double tol)
{
    typedef typename VECT::SCALAR S;
    typedef typename alg_types<2,2,Rational>::S RatType;
    check_vec_almost_equal_impl<boost::is_same<S, RatType>::value>::
        check(expected, actual, tol);
}


#define CHECK_VEC_CLOSE(expected, actual, tol) check_vec_almost_equal(expected, actual, tol)


#endif