//
// Created by sam on 26/10/2021.
//

#include <UnitTest++/UnitTest++.h>
#include <libalgebra/libalgebra.h>
#include <libalgebra/coefficients/coefficients.h>

#include "fixture.h"

struct Fixture : fixture_base {


    using basis_t = alg::free_tensor_basis<5, 5>;
    using coeff_t = alg::coefficients::double_field;

    using vector_t = alg::vectors::dense_vector<basis_t, coeff_t>;

};


SUITE(dense_double_serialization) {

    TEST_FIXTURE(Fixture, test_serialize_empty_vec)
    {
        vector_t vec;


    }


}