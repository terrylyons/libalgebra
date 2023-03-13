//
// Created by sam on 28/10/2021.
//

#include <UnitTest++.h>
#include <libalgebra/alg_types.h>
#include <libalgebra/coefficients.h>
#include <libalgebra/rational_coefficients.h>
#include <libalgebra/libalgebra.h>

#include "../../common/random_vector_generator.h"
#include "../../common/rng.h"
#include "fixture.h"

SUITE(higher_class_double_serialization)
{

    struct Fixture : fixture_base, alg_types<5, 5, DPReal> {
        using coeff_t = alg_types<5, 5, DPReal>::COEFF;
        using skip_dist = std::uniform_int_distribution<alg::DIMN>;
        using coeff_dist = std::uniform_real_distribution<typename coeff_t::S>;

        skip_dist sd;
        coeff_dist cd;

        template<typename Vector>
        using rvg_t = la_testing::random_vector_generator<Vector, coeff_dist>;
        std::mt19937 rng;

        Fixture() : fixture_base(), sd(0, 5), cd(-1.0, 1.0), rng(std::random_device()())
        {}
    };

    TEST_FIXTURE(Fixture, serialize_free_tensor)
    {
        auto ft = rvg_t<TENSOR>(-1.0, 1.0)(rng);

        write(ft);
        auto read_vec = read<TENSOR>();

        CHECK_EQUAL(ft, read_vec);
    }

    TEST_FIXTURE(Fixture, serialize_shuffle_tensor)
    {
        auto ft = rvg_t<SHUFFLE_TENSOR>(-1.0, 1.0)(rng);

        write(ft);
        auto read_vec = read<SHUFFLE_TENSOR>();

        CHECK_EQUAL(ft, read_vec);
    }

    TEST_FIXTURE(Fixture, serialize_lie)
    {
        auto ft = rvg_t<LIE>(-1.0, 1.0)(rng);

        write(ft);
        auto read_vec = read<LIE>();

        CHECK_EQUAL(ft, read_vec);
    }
}

SUITE(higher_class_rational_serialization)
{

    struct Fixture : fixture_base, alg_types<5, 5, Rational> {
        using coeff_t = alg_types<5, 5, Rational>::COEFF;
        using skip_dist = std::uniform_int_distribution<alg::DIMN>;
        using coeff_dist = la_testing::uniform_rational_distribution<typename coeff_t::S>;

        skip_dist sd;
        coeff_dist cd;

        template<typename Vector>
        using rvg_t = la_testing::random_vector_generator<Vector, coeff_dist>;
        std::mt19937 rng;

        Fixture() : fixture_base(), sd(0, 5), cd(-1.0, 1.0), rng(std::random_device()())
        {}
    };

    TEST_FIXTURE(Fixture, serialize_free_tensor)
    {
        auto ft = rvg_t<TENSOR>(S(-1), S(1))(rng);

        write(ft);
        auto read_vec = read<TENSOR>();

        CHECK_EQUAL(ft, read_vec);
    }

    TEST_FIXTURE(Fixture, serialize_shuffle_tensor)
    {
        auto ft = rvg_t<SHUFFLE_TENSOR>(S(-1), S(1))(rng);

        write(ft);
        auto read_vec = read<SHUFFLE_TENSOR>();

        CHECK_EQUAL(ft, read_vec);
    }

    TEST_FIXTURE(Fixture, serialize_lie)
    {
        auto ft = rvg_t<LIE>(S(-1), S(1))(rng);

        write(ft);
        auto read_vec = read<LIE>();

        CHECK_EQUAL(ft, read_vec);
    }
}
