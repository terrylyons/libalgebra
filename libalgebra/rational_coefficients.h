//
// Created by sam on 15/10/2021.
//

#ifndef LIBALGEBRA_RATIONAL_COEFFICIENTS_H
#define LIBALGEBRA_RATIONAL_COEFFICIENTS_H

#ifndef LIBALGEBRA_NO_GMP
#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
#include "coefficients/gmp_ser.h"

#endif
#include <boost/multiprecision/gmp.hpp>
#else
#include <boost/multiprecision/cpp_int.hpp>
#endif

#include "coefficients.h"

namespace alg {
namespace coefficients {

namespace mp = boost::multiprecision;

#ifndef LIBALGEBRA_NO_GMP
typedef mp::gmp_rational rational_backend;
#else
typedef mp::cpp_rational_backend rational_backend;
#endif

typedef mp::number<rational_backend, mp::et_off> rational;

typedef coefficient_field<rational> rational_field;

}// namespace coefficients
}// namespace alg

#endif//LIBALGEBRA_RATIONAL_COEFFICIENTS_H
