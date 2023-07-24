/* *************************************************************

Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurkó and Arend Janssen.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */


//  implimentation_types.h : provides definitions for basic types

#ifndef implimentation_types_h__
#define implimentation_types_h__

#include <stddef.h>
#include <type_traits>

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

//#define ORDEREDMAP
#define UNORDEREDMAP
#define NOBTREE

#if (!(defined _MSC_VER) && (__cplusplus < 201103L)) || ((defined _MSC_VER) && (_MSC_VER < 1800))
// we do not not have a C++11 compliant compiler
// visual studio 2008 does not compile the btree or have unordered map header or
// have variadic templates gcc 4.8 does not support C++11 unless it is switched
// on
//#if (_MSC_VER < 1800) // VS2008 VER 1500 is still needed for python 2.7 VS2010
//VER 1700 is still needed for python 3.4
#ifndef NOBTREE
#define NOBTREE
#endif// !NOBTREE
#ifndef ORDEREDMAP
#define ORDEREDMAP
#endif// !ORDEREDMAP
#ifdef UNORDEREDMAP
#undef UNORDEREDMAP
#endif// UNORDEREDMAP
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
#endif// !UNORDEREDMAP
#ifndef NOBTREE
#include "cpp-btree/safe_btree_map.h"
#endif// !NOBTREE


namespace alg {

// moved from  libalgebra.h

/// Used to store degrees.
typedef unsigned DEG;
typedef typename std::make_signed<DEG>::type IDEG;// signed integer same size as degree

/// Used to index letters, and basis elements. The value 0 may be special.
typedef size_t LET;

/// Used for large integer calculations where overflow might otherwise occur
typedef unsigned long long LET64;

/// Used for dimension of of vector elements and indices in dense vectors
typedef size_t DIMN;
typedef std::ptrdiff_t IDIMN;

}// namespace alg
#endif// implimetation_types_h__
