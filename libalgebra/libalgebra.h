/* *************************************************************

Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurk� and Arend Janssen.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */

//  libalgebra.h

// Include once wrapper
#ifndef DJC_COROPA_LIBALGEBRAH_SEEN
#define DJC_COROPA_LIBALGEBRAH_SEEN

#ifdef _MSC_VER
#if _MSC_VER >= 1600
#include <cstdint>
#else
typedef __int8 int8_t;
typedef __int16 int16_t;
typedef __int32 int32_t;
typedef __int64 int64_t;
typedef unsigned __int8 uint8_t;
typedef unsigned __int16 uint16_t;
typedef unsigned __int32 uint32_t;
typedef unsigned __int64 uint64_t;
#endif
#elif __cplusplus < 201103L

#include <stdint.h>

#else

#include <cstdint>

#endif

#include <algorithm>
#include <cassert>
#include <cmath>
#include <deque>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <utility>
#include <vector>

#include <boost/thread/recursive_mutex.hpp>
#include <boost/thread/shared_mutex.hpp>

#include "compat.h"
#include "implimentation_types.h"
#include <boost/thread/locks.hpp>

#if __cplusplus >= 201103L

#include <array>

#define LIBALGEBRA_STATIC_ARRAY_TYPE std::array
#else
#include <boost/array.hpp>
#define LIBALGEBRA_STATIC_ARRAY_TYPE boost::array
#endif

//#define ORDEREDMAP
#define UNORDEREDMAP
#define NOBTREE

#if (!(defined _MSC_VER) && (__cplusplus < 201103L)) || \
    ((defined _MSC_VER) && (_MSC_VER < 1800))
// we do not not have a C++11 compliant compiler
// visual studio 2008 does not compile the btree or have unordered map header or
// have variadic templates gcc 4.8 does not support C++11 unless it is switched
// on
//#if (_MSC_VER < 1800) // VS2008 VER 1500 is still needed for python 2.7 VS2010
//VER 1700 is still needed for python 3.4
#ifndef NOBTREE
#define NOBTREE
#endif // !NOBTREE
#ifndef ORDEREDMAP
#define ORDEREDMAP
#endif // !ORDEREDMAP
#ifdef UNORDEREDMAP
#undef UNORDEREDMAP
#endif // UNORDEREDMAP
//#endif
//#endif
#endif

#ifdef UNORDEREDMAP
// require C++11 support
//#include "addons/sized_unordered_map.h"
//#define MY_UNORDERED_MAP sized_unordered_map
#include <type_traits>
#include <unordered_map>

#define MY_UNORDERED_MAP std::unordered_map
#else
#define ORDEREDMAP
#endif // !UNORDEREDMAP
#ifndef NOBTREE
#include "cpp-btree/safe_btree_map.h"
#endif // !NOBTREE

/// The libalgebra namespace. A set of template classes for Algebras.
/**
   The SCA typename corresponds to a commutative ring with unit containing the
   rationals. Needed operators are * + - and explicit ctor from int type.

   The RAT typename corresponds to rationals built from SCA. Needed
   operators are + - * / and explicit ctor from SCA type and int type.
*/

#include <boost/container/flat_map.hpp>

#include "libalgebra/utils/integer_maths.h"
#include "libalgebra/vectors/vectors.h"
#include "multiplication_helpers.h"



namespace alg {
#include "constlog2.h"
#include "constpower.h"
}
// Some useful macros to avoid similar codes.

#define __DECLARE_BINARY_OPERATOR(T1, NEWOP, OLDOP, T2)                        \
  T1 operator NEWOP(const T2 &rhs) const {                                     \
    T1 result(*this);                                                          \
    result OLDOP rhs;                                                          \
    return result;                                                             \
  }

#define __DECLARE_UNARY_OPERATOR(NEWT, NEWOP, OLDOP, OLDT)                     \
  NEWT operator NEWOP(void) const { return OLDT::operator OLDOP(); }

//  End of macros.

/// Forward declaration of classes

/// Sparse vectors with default MAP typename from BASIS typename.
// template<class BASIS, class MAP = typename BASIS::MAP>
// class sparse_vector;

namespace alg {
#include "algebra.h"

#include "lie_basis.h"
#include "monomial_basis.h"
#include "poly_basis.h"
#include "poly_lie_basis.h"
}

#include "tensor_basis.h"


namespace alg {
#include "lie.h"
#include "multi_polynomial.h"
#include "poly_lie.h"
#include "polynomials.h"
#include "tensor.h"



#include "utils.h"
}
// Undeclaring local macros in reverse order of declaration.
#undef __DECLARE_UNARY_OPERATOR
#undef __DECLARE_BINARY_OPERATOR
// End of undeclaration of local macros.



// Include once wrapper
#endif // DJC_COROPA_LIBALGEBRAH_SEEN

// EOF.
