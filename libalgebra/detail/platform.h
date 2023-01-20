//
// Created by user on 19/01/23.
//

#ifndef LIBALGEBRA_LIBALGEBRA_DETAIL_PLATFORM_H_
#define LIBALGEBRA_LIBALGEBRA_DETAIL_PLATFORM_H_

#include <boost/predef.h>

//#define LIBALGEBRA_ALLOW_PREFETCHING

#ifdef LIBALGEBRA_ALLOW_PREFETCHING


#if defined(__has_include) && __has_include("xmmintrin.h")
#include <xmmintrin.h>

#define LA_PREFETCH_T0(PTR) _mm_prefetch((PTR), _MM_HINT_T0)
#define LA_PREFETCH_T1(PTR) _mm_prefetch((PTR), _MM_HINT_T1)
#define LA_PREFETCH_T2(PTR) _mm_prefetch((PTR), _MM_HINT_T2)
#ifdef _MM_HINT_ET0
#define LA_PREFETCH_ET0(PTR) _mm_prefetch((PTR), _MM_HINT_ET0)
#else
#define LA_PREFETCH_ET0(PTR) _mm_prefetch((PTR), _MM_HINT_T0)
#endif
#ifdef _MM_HINT_ET1
#define LA_PREFETCH_ET1(PTR) _mm_prefetch((PTR), _MM_HINT_ET1)
#else
#define LA_PREFETCH_ET1(PTR) _mm_prefetch((PTR), _MM_HINT_T1)
#endif


#endif

#endif

#ifndef LA_PREFETCH_T0
#define LA_PREFETCH_T0(PTR)
#define LA_PREFETCH_T1(PTR)
#define LA_PREFETCH_T2(PTR)
#define LA_PREFETCH_ET0(PTR)
#define LA_PREFETCH_ET1(PTR)
#endif



#endif//LIBALGEBRA_LIBALGEBRA_DETAIL_PLATFORM_H_
