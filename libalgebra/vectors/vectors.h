//
// Created by sam on 01/02/2021.
//

#ifndef LIBALGEBRA_VECTORS_H
#define LIBALGEBRA_VECTORS_H


#define __DECLARE_BINARY_OPERATOR(T1, NEWOP, OLDOP, T2) \
	T1 operator NEWOP(const T2& rhs) const \
	{ T1 result(*this); return result OLDOP rhs; }

#define __DECLARE_UNARY_OPERATOR(NEWT, NEWOP, OLDOP, OLDT) \
	NEWT operator NEWOP(void) const \
	{ return OLDT::operator OLDOP (); }

#include "libalgebra/vectors/sparse_vector.h"
#include "libalgebra/vectors/dense_vector.h"

#include "libalgebra/vectors/vector.h"

#undef __DECLARE_UNARY_OPERATOR
#undef __DECLARE_BINARY_OPERATOR


#endif //LIBALGEBRA_VECTORS_H
