//
// Created by sam on 02/02/2021.
//


#include <UnitTest++.h>

#include <libalgebra/libalgebra.h>
#include <libalgebra/alg_types.h>

#include "../../common/time_and_details.h"

#define CHECK_VEC_CLOSE(expected, result, tol) \
    CHECK_EQUAL(expected, result)

using alg::DEG;
using alg::LET;

struct Fixture : alg_types<5, 5, DPReal>
{
    typedef alg_types<5, 5, DPReal> AT;
    typedef typename AT::S S;
    typedef typename AT::Q Q;
    typedef typename AT::LIE LIE;
    typedef typename LIE::KEY TKEY;
};


SUITE(lie_basis) {





}
