
/* *************************************************************

Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai, 
Greg Gyurk� and Arend Janssen. 

Distributed under the terms of the GNU General Public License, 
Version 3. (See accompanying file License.txt)

************************************************************* */


#include <stddef.h>

//  implimentation_types.h : provides definitions for basic types

#ifndef implimentation_types_h__
#define implimentation_types_h__

namespace alg {

// moved from  libalgebra.h

/// Used to store degrees.
typedef unsigned DEG;
typedef int IDEG; // signed integer same size as degree

/// Used to index letters, and basis elements. The value 0 may be special.
typedef size_t LET;

/// Used for large integer calculations where overflow might otherwise occur
typedef unsigned long long LET64;

/// Used for dimension of of vector elements and indices in dense vectors
typedef size_t DIMN;
typedef std::ptrdiff_t IDIMN;

}
#endif // implimetation_types_h__