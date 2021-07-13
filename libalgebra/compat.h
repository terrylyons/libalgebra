//
// Created by sam on 13/07/2021.
//

#ifndef LIBALGEBRA_COMPAT_H
#define LIBALGEBRA_COMPAT_H


#if __cplusplus >= 201103UL
#define LA_CONSTEXPR constexpr
#define LA_EXPLICIT explicit
#else
#define LA_CONSTEXPR
#define LA_EXPLICIT
#endif






#endif //LIBALGEBRA_COMPAT_H
