//
// Created by sam on 26/10/2021.
//

#include <UnitTest++/UnitTest++.h>
#include <libalgebra/coefficients/coefficients.h>
#include <libalgebra/coefficients/rational_coefficients.h>
#include <libalgebra/libalgebra.h>

#include "../../common/random_vector_generator.h"
#include "fixture.h"

SUITE(dense_double_serialization)
{
    struct Fixture : fixture_base {

        using basis_t = alg::free_tensor_basis<5, 3>;
        using coeff_t = alg::coefficients::double_field;
        using key_type = typename basis_t::KEY;

        using vector_t = alg::vectors::dense_vector<basis_t, coeff_t>;
    };

    TEST_FIXTURE(Fixture, test_serialize_empty_vec)
    {
        vector_t vec;

        write(vec);
        auto read_vec = read<vector_t>();

        CHECK_EQUAL(vec, read_vec);
    }

    TEST_FIXTURE(Fixture, test_serialize_random_vec)
    {
        la_testing::random_vector_generator<vector_t> rvg(-1.0, 1.0);
        boost::random::mt19937 rng;

        vector_t vec = rvg(rng);
        REQUIRE CHECK_EQUAL(basis_t::start_of_degree(4), vec.size());

        write(vec);
        auto read_vec = read<vector_t>();

        CHECK_EQUAL(vec, read_vec);
    }
}

SUITE(dense_float_serialization)
{
    struct Fixture : fixture_base {

        using basis_t = alg::free_tensor_basis<5, 3>;
        using coeff_t = alg::coefficients::float_field;
        using key_type = typename basis_t::KEY;

        using vector_t = alg::vectors::dense_vector<basis_t, coeff_t>;
    };

    TEST_FIXTURE(Fixture, test_serialize_empty_vec)
    {
        vector_t vec;

        write(vec);
        auto read_vec = read<vector_t>();

        CHECK_EQUAL(vec, read_vec);
    }

    TEST_FIXTURE(Fixture, test_serialize_random_vec)
    {
        la_testing::random_vector_generator<vector_t> rvg(-1.0f, 1.0f);
        boost::random::mt19937 rng;

        vector_t vec = rvg(rng);
        REQUIRE CHECK_EQUAL(basis_t::start_of_degree(4), vec.size());

        write(vec);
        auto read_vec = read<vector_t>();

        CHECK_EQUAL(vec, read_vec);
    }
}

/*
SUITE(dense_rational_serialization)
{
    struct Fixture : fixture_base {

        using basis_t = alg::free_tensor_basis<5, 3>;
        using coeff_t = alg::coefficients::rational_field;
        using key_type = typename basis_t::KEY;

        using vector_t = alg::vectors::dense_vector<basis_t, coeff_t>;
    };

    TEST_FIXTURE(Fixture, test_serialize_empty_vec)
    {
        vector_t vec;

        write(vec);
        auto read_vec = read<vector_t>();

        CHECK_EQUAL(vec, read_vec);
    }

    TEST_FIXTURE(Fixture, test_serialize_random_vec)
    {
        la_testing::random_vector_generator<vector_t> rvg(-1.0, 1.0);
        boost::random::mt19937 rng;

        vector_t vec = rvg(rng);
        REQUIRE CHECK_EQUAL(basis_t::start_of_degree(4), vec.size());

        write(vec);
        auto read_vec = read<vector_t>();

        CHECK_EQUAL(vec, read_vec);
    }
}
*/