//
// Created by sam on 26/10/2021.
//
#include <UnitTest++.h>
#include <libalgebra/coefficients.h>
#include <libalgebra/rational_coefficients.h>
#include <libalgebra/libalgebra.h>

#include <map>
#include <unordered_map>

#include "../../common/random_vector_generator.h"
#include "../../common/rng.h"
#include "fixture.h"

SUITE(sparse_map_double_serialization)
{
    struct Fixture : fixture_base {

        using basis_t = alg::free_tensor_basis<5, 3>;
        using coeff_t = alg::coefficients::double_field;
        using key_type = typename basis_t::KEY;
        using map_t = std::map<key_type, typename coeff_t::S>;

        using vector_t = alg::vectors::sparse_vector<basis_t, coeff_t, map_t>;

        using skip_dist = std::uniform_int_distribution<alg::DIMN>;
        using coeff_dist = std::uniform_real_distribution<typename coeff_t::S>;

        skip_dist sd;
        coeff_dist cd;

        using rvg_t = la_testing::random_vector_generator<vector_t, coeff_dist, skip_dist>;
        rvg_t rvg;
        std::mt19937 rng;

        Fixture() : fixture_base(), sd(0, 5), cd(-1.0, 1.0),
                    rvg(sd, -1.0, 1.0), rng(std::random_device()())
        {}
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
        vector_t vec = rvg(rng);

        write(vec);
        auto read_vec = read<vector_t>();

        CHECK_EQUAL(vec, read_vec);
    }
}

SUITE(sparse_map_float_serialization)
{
    struct Fixture : fixture_base {

        using basis_t = alg::free_tensor_basis<5, 3>;
        using coeff_t = alg::coefficients::float_field;
        using key_type = typename basis_t::KEY;
        using map_t = std::map<key_type, typename coeff_t::S>;

        using vector_t = alg::vectors::sparse_vector<basis_t, coeff_t, map_t>;

        using skip_dist = std::uniform_int_distribution<alg::DIMN>;
        using coeff_dist = std::uniform_real_distribution<typename coeff_t::S>;

        skip_dist sd;
        coeff_dist cd;

        using rvg_t = la_testing::random_vector_generator<vector_t, coeff_dist, skip_dist>;
        rvg_t rvg;
        std::mt19937 rng;

        Fixture() : fixture_base(), sd(0, 5), cd(-1.0f, 1.0f),
                    rvg(sd, -1.0f, 1.0f), rng(std::random_device()())
        {}
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
        vector_t vec = rvg(rng);

        write(vec);
        auto read_vec = read<vector_t>();

        CHECK_EQUAL(vec, read_vec);
    }
}

SUITE(sparse_map_rational_serialization)
{
    struct Fixture : fixture_base {

        using basis_t = alg::free_tensor_basis<5, 3>;
        using coeff_t = alg::coefficients::rational_field;
        using key_type = typename basis_t::KEY;
        using map_t = std::map<key_type, typename coeff_t::S>;

        using vector_t = alg::vectors::sparse_vector<basis_t, coeff_t, map_t>;

        using skip_dist = std::uniform_int_distribution<alg::DIMN>;
        using coeff_dist = la_testing::uniform_rational_distribution<typename coeff_t::S>;

        skip_dist sd;
        coeff_dist cd;

        using rvg_t = la_testing::random_vector_generator<vector_t, coeff_dist, skip_dist>;
        rvg_t rvg;
        std::mt19937 rng;

        Fixture() : fixture_base(), sd(0, 5), cd(-1, 1),
                    rvg(sd, -1.0, 1.0), rng(std::random_device()())
        {}
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
        vector_t vec = rvg(rng);

        write(vec);
        auto read_vec = read<vector_t>();

        CHECK_EQUAL(vec, read_vec);
    }
}

SUITE(sparse_unordered_map_double_serialization)
{
    struct Fixture : fixture_base {

        using basis_t = alg::free_tensor_basis<5, 3>;
        using coeff_t = alg::coefficients::double_field;
        using key_type = typename basis_t::KEY;
        using map_t = std::unordered_map<key_type, typename coeff_t::S>;

        using vector_t = alg::vectors::sparse_vector<basis_t, coeff_t, map_t>;

        using skip_dist = std::uniform_int_distribution<alg::DIMN>;
        using coeff_dist = std::uniform_real_distribution<typename coeff_t::S>;

        skip_dist sd;
        coeff_dist cd;

        using rvg_t = la_testing::random_vector_generator<vector_t, coeff_dist, skip_dist>;
        rvg_t rvg;
        std::mt19937 rng;

        Fixture() : fixture_base(), sd(0, 5), cd(-1.0, 1.0),
                    rvg(sd, -1.0, 1.0), rng(std::random_device()())
        {}
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
        vector_t vec = rvg(rng);

        write(vec);
        auto read_vec = read<vector_t>();

        CHECK_EQUAL(vec, read_vec);
    }
}

SUITE(sparse_unordered_map_float_serialization)
{
    struct Fixture : fixture_base {

        using basis_t = alg::free_tensor_basis<5, 3>;
        using coeff_t = alg::coefficients::float_field;
        using key_type = typename basis_t::KEY;
        using map_t = std::unordered_map<key_type, typename coeff_t::S>;

        using vector_t = alg::vectors::sparse_vector<basis_t, coeff_t, map_t>;

        using skip_dist = std::uniform_int_distribution<alg::DIMN>;
        using coeff_dist = std::uniform_real_distribution<typename coeff_t::S>;

        skip_dist sd;
        coeff_dist cd;

        using rvg_t = la_testing::random_vector_generator<vector_t, coeff_dist, skip_dist>;
        rvg_t rvg;
        std::mt19937 rng;

        Fixture() : fixture_base(), sd(0, 5), cd(-1.0f, 1.0f),
                    rvg(sd, -1.0f, 1.0f), rng(std::random_device()())
        {}
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
        vector_t vec = rvg(rng);

        write(vec);
        auto read_vec = read<vector_t>();

        CHECK_EQUAL(vec, read_vec);
    }
}

SUITE(sparse_unordered_map_rational_serialization)
{
    struct Fixture : fixture_base {

        using basis_t = alg::free_tensor_basis<5, 3>;
        using coeff_t = alg::coefficients::rational_field;
        using key_type = typename basis_t::KEY;
        using map_t = std::unordered_map<key_type, typename coeff_t::S>;

        using vector_t = alg::vectors::sparse_vector<basis_t, coeff_t, map_t>;

        using skip_dist = std::uniform_int_distribution<alg::DIMN>;
        using coeff_dist = la_testing::uniform_rational_distribution<typename coeff_t::S>;

        skip_dist sd;
        coeff_dist cd;

        using rvg_t = la_testing::random_vector_generator<vector_t, coeff_dist, skip_dist>;
        rvg_t rvg;
        std::mt19937 rng;

        Fixture() : fixture_base(), sd(0, 5), cd(-1, 1),
                    rvg(sd, -1.0, 1.0), rng(std::random_device()())
        {}
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
        vector_t vec = rvg(rng);

        write(vec);
        auto read_vec = read<vector_t>();

        CHECK_EQUAL(vec, read_vec);
    }
}
