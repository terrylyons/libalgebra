//
// Created by sam on 13/07/2021.
//

#ifndef LIBALGEBRA_COMPAT_H
#define LIBALGEBRA_COMPAT_H


#if __cplusplus >= 201103UL
#define LA_CONSTEXPR constexpr
#else
#define LA_CONSTEXPR
#endif






#endif //LIBALGEBRA_COMPAT_H
