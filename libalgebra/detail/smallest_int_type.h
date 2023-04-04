//
// Created by user on 05/12/22.
//

#ifndef LIBALGEBRA_LIBALGEBRA_DETAIL_SMALLEST_INT_TYPE_H_
#define LIBALGEBRA_LIBALGEBRA_DETAIL_SMALLEST_INT_TYPE_H_


#include "integer_maths.h"
#include <boost/integer.hpp>


namespace alg {
namespace dtl {

template <int Digits>
struct smallest_containing_int
{
    static constexpr int num_bits = ::alg::integer_maths::logN(Digits, 2) + 1;
    using type = typename boost::uint_t<num_bits>::least;

};




} // namespace dtl
} // namespace alg

#endif//LIBALGEBRA_LIBALGEBRA_DETAIL_SMALLEST_INT_TYPE_H_
