//
// Created by sam on 03/02/2021.
//


#include <map>

#include <UnitTest++.h>

#include "../../common/simple_basis.h"
#include "../../common/time_and_details.h"
#include <libalgebra/alg_types.h>

#include "../../common/compat.h"
#include "../../common/rng.h"


#ifdef LIBALGEBRA_VECTORS_H
using alg::vectors::sparse_vector;
#else
using alg::sparse_vector;
#endif

#define _VECTOR_TYPE sparse_vector
#define _VECTOR_TYPE_ADDITIONAL_PARAMS std::map<typename BASIS::KEY, typename FIELD::S>
#define _TVECTOR_TYPE_ADDITIONAL_PARAMS \
    std::map<typename TBASIS::KEY, typename FIELD::S>

SUITE(sparse_vector) {

#include "framework_fixture.h"
#include "vector_arithmetic_suite.h"
#include "vector_comparison_suite.h"
#include "vector_element_access_suite.h"
#include "vector_iterator_suite.h"
#include "vector_properties_suite.h"
#include "vector_norm_suite.h"
//#include "vector_binary_transform_suite.h"

}


#undef _VECTOR_TYPE_ADDITIONAL_PARAMS
#undef _VECTOR_TYPE
