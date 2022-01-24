//
// Created by sam on 17/11/2021.
//

#include "basis_iteration_suite.h"
#include <libalgebra/libalgebra.h>

static auto tensor22_suite = create_basis_test_suite<alg::tensor_basis<2, 2>>("tensor_basis_22");
static auto tensor35_suite = create_basis_test_suite<alg::tensor_basis<2, 5>>("tensor_basis_35");
static auto tensor102_suite = create_basis_test_suite<alg::tensor_basis<10, 2>>("tensor_basis_102");
