//
// Created by sam on 26/10/2021.
//

#include <UnitTest++/UnitTest++.h>
#include <libalgebra/coefficients/rational_coefficients.h>
#include <libalgebra/vectors/dense_storage.h>
#include <random>

#include "fixture.h"

struct fixture : fixture_base {
    using storage = alg::vectors::dense_storage<double>;
    using size_type = storage::size_type;
};

SUITE(dense_storage_serialization)
{

    TEST_FIXTURE(fixture, test_empty_storage_serialize)
    {
        storage x;

        write(x);

        storage y = read<storage>();

        CHECK_EQUAL(size_type(0), y.size());
        CHECK_EQUAL(x, y);
    }

    TEST_FIXTURE(fixture, test_non_empty_storage_serialize)
    {
        storage x{1.0, 2.0, 3.0, 4.0, 5.0};

        write(x);

        storage y = read<storage>();

        CHECK_EQUAL(size_type(5), y.size());
        CHECK_EQUAL(x, y);
    }

    TEST_FIXTURE(fixture, test_serialization_borrowed_random)
    {
        std::vector<double> data;
        data.reserve(50);

        std::random_device rd;
        std::mt19937 rng(rd());
        std::uniform_real_distribution<double> dist(-1.0, 1.0);

        for (int i = 0; i < 50; ++i) {
            data.emplace_back(dist(rng));
        }

        const double* begin = data.data();
        const double* end = begin + data.size();
        storage borrowed(begin, end);
        REQUIRE CHECK_EQUAL(size_type(50), borrowed.size());
        REQUIRE CHECK_EQUAL(storage::vec_type::borrowed, borrowed.type());

        write(borrowed);

        storage y = read<storage>();

        CHECK_EQUAL(size_type(50), y.size());
        CHECK_EQUAL(storage::vec_type::owned, y.type());
        CHECK_EQUAL(borrowed, y);
    }

    TEST_FIXTURE(fixture, test_serialization_borrowed_mut_random)
    {
        std::vector<double> data;
        data.reserve(50);

        std::random_device rd;
        std::mt19937 rng(rd());
        std::uniform_real_distribution<double> dist(-1.0, 1.0);

        for (int i = 0; i < 50; ++i) {
            data.emplace_back(dist(rng));
        }
        
        double* begin = data.data();
        double* end = begin + data.size();
        storage borrowed(begin, end);
        REQUIRE CHECK_EQUAL(size_type(50), borrowed.size());
        REQUIRE CHECK_EQUAL(storage::vec_type::borrowed_mut, borrowed.type());

        write(borrowed);

        storage y = read<storage>();

        CHECK_EQUAL(size_type(50), y.size());
        CHECK_EQUAL(storage::vec_type::owned, y.type());
        CHECK_EQUAL(borrowed, y);
    }

    TEST_FIXTURE(fixture, test_serialization_rational_storage)
    {
        using rational = alg::coefficients::rational;
        alg::vectors::dense_storage<rational> data{rational(1), rational(2), rational(3.1415)};

        write(data);

        auto y = read<alg::vectors::dense_storage<rational>>();

        CHECK_EQUAL(data, y);
    }
}