//
// Created by sam on 23/02/2022.
//

#include <UnitTest++/UnitTest++.h>
#include <libalgebra/alg_types.h>
#include <libalgebra/vectors/vector.h>

#include <vector>

SUITE(tensor_creation) {

    struct fixture : alg_types<5, 5, DPReal>
    {
        using dense_tensor = alg::free_tensor<COEFF, 5, 5, alg::vectors::dense_vector>;
    };

    TEST_FIXTURE(fixture, test_creation_borrowed_data_resize) {
        std::vector<double> values {1.0, 2.0}; // Definitely not a full dimension
        dense_tensor constructed(values.data(), values.data()+2);

        CHECK_EQUAL(6, alg::vectors::dtl::vector_base_access::convert(constructed).dimension());
    }

}