//
// Created by sam on 03/02/2021.
//
#include "vector_suite.h"

#include <UnitTest++/UnitTest++.h>

#include "../../common/compat.h"
#include "../../common/rng.h"
#include "../../common/simple_basis.h"
#include "../../common/time_and_details.h"
#include <libalgebra/alg_types.h>
#include <libalgebra/libalgebra.h>

#ifdef LIBALGEBRA_VECTORS_H
using alg::vectors::dense_vector;

#define _VECTOR_TYPE dense_vector
#undef _VECTOR_TYPE_ADDITIONAL_PARAMS

SUITE(dense_vector)
{

#include "framework_fixture.h"
#include "vector_arithmetic_suite.h"
#include "vector_binary_transform_suite.h"
#include "vector_comparison_suite.h"
#include "vector_element_access_suite.h"
#include "vector_iterator_suite.h"
#include "vector_norm_suite.h"
#include "vector_properties_suite.h"
}

#undef _VECTOR_TYPE_ADDITIONAL_PARAMS
#undef _VECTOR_TYPE
#endif

namespace {
MAKE_SUITE_FOR(test_dense_vector, vector_suite, alg::vectors::dense_vector);
}