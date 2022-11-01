//
// Created by sam on 12/02/2021.
//

#include <UnitTest++.h>

#include <libalgebra/libalgebra.h>
#include <libalgebra/alg_types.h>
#include "../../common/simple_basis.h"
#include "../../common/time_and_details.h"
#include "../../common/rng.h"

using alg::DEG;


SUITE(dense_algebra) {
#define ALGEBRA_TESTS_VECT_TYPE alg::vectors::dense_vector

#include "simple_algebra_basis.h"
#include "framework_fixture.h"

#include "alg_mul_suite.h"


}
#undef ALGEBRA_TESTS_VECT_TYPE

SUITE(sparse_algebra) {
#define ALGEBRA_TESTS_VECT_TYPE \
    alg::vectors::sparse_vector

#include "simple_algebra_basis.h"
#include "framework_fixture.h"

#include "alg_mul_suite.h"


}
#undef ALGEBRA_TESTS_VECT_TYPE
