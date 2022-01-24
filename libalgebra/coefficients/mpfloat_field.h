//
// Created by sam on 24/01/2022.
//

#ifndef LIBALGEBRA_MPFLOAT_FIELD_H
#define LIBALGEBRA_MPFLOAT_FIELD_H

#include <libalgebra/implementation_types.h>

#ifndef LIBALGEBRA_NO_GMP
#include <boost/multiprecision/gmp.hpp>
#include <libalgebra/coefficients/gmp_ser.h>
#else
#include <boost/multipricision/cpp_bin_float.hpp>
#endif

#include <libalgebra/coefficients/coefficients.h>

namespace alg {
namespace coefficients {

namespace mp = boost::multiprecision;

#ifndef LIBALGEBRA_NO_GMP
template<DEG Digits10>
using mpfloat_backend = mp::gmp_float<Digits10>;
#else
template<DEG Digits10>
using mpfloat_backend = mp::cpp_bin_float<Digits10>;
#endif

template<DEG Digits10>
using mpfloat = mp::number<mpfloat_backend<Digits10>, mp::et_off>;

template<DEG Digits10>
using mpfloat_field = coefficient_field<mpfloat<Digits10>>;

}// namespace coefficients
}// namespace alg

#endif//LIBALGEBRA_MPFLOAT_FIELD_H
