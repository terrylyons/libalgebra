//
// Created by sam on 26/10/2021.
//

#include <UnitTest++/UnitTest++.h>
#include <libalgebra/vectors/dense_storage.h>

#include "fixture.h"

struct fixture : fixture_base {
    using storage = alg::vectors::dense_storage<double>;
};

SUITE(dense_storage_serialization)
{

    TEST_FIXTURE(fixture, test_empty_storage_serialize)
    {
        storage x;

        write(x);

        storage y = read<storage>();

        CHECK_EQUAL(x, y);
    }

    TEST_FIXTURE(fixture, test_non_empty_storage_serialize)
    {
        storage x{1.0, 2.0, 3.0, 4.0, 5.0};

        write(x);

        storage y = read<storage>();

        CHECK_EQUAL(x, y);
    }






}