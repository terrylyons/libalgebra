//
// Created by sam on 17/11/2021.
//

#include "basis_iteration_suite.h"
#include <libalgebra/libalgebra.h>

static auto lie22_suite = create_basis_test_suite<alg::lie_basis<2, 2>>("lie_basis_22");
static auto lie35_suite = create_basis_test_suite<alg::lie_basis<2, 5>>("lie_basis_35");
static auto lie102_suite = create_basis_test_suite<alg::lie_basis<10, 2>>("lie_basis_102");
