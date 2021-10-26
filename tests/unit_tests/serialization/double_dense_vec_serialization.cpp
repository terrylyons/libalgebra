//
// Created by sam on 26/10/2021.
//

#include <UnitTest++/UnitTest++.h>
#include <libalgebra/coefficients/coefficients.h>
#include <libalgebra/libalgebra.h>

#include "fixture.h"

struct Fixture : fixture_base {

    using basis_t = alg::free_tensor_basis<5, 5>;
    using coeff_t = alg::coefficients::double_field;
    using key_type = typename basis_t::KEY;

    using vector_t = alg::vectors::dense_vector<basis_t, coeff_t>;
};

SUITE(dense_double_serialization)
{

    TEST_FIXTURE(Fixture, test_serialize_empty_vec)
    {
        vector_t vec;

        write(vec);
        auto read_vec = read<vector_t>();

        CHECK_EQUAL(vec, read_vec);
    }
}