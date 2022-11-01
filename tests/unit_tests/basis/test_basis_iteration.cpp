//
// Created by sam on 16/11/2021.
//

#include <UnitTest++.h>
#include <libalgebra/key_iterators.h>
#include <libalgebra/libalgebra.h>


#ifdef __GNUC__
#define UNUSED __attribute__((unused))
#elif defined(_MSC_VER)
#define UNUSED __pragma(warning(suppress:4100))
#else
#define UNUSED
#endif

SUITE(test_basis_iterator_helpers)
{

    struct fixture {
        using tensor_basis = alg::tensor_basis<2, 2>;
        using key_type = typename tensor_basis::KEY;
        using key_range = alg::basis::key_range<tensor_basis>;
        using iterator = alg::basis::key_iterator<tensor_basis>;

        tensor_basis basis;

        key_range iterate_basis() const
        {
            return key_range(basis);
        }

        key_range iterate_basis(const key_type& k) const
        {
            return key_range(basis, k);
        }

        key_range iterate_basis(const key_type& start, const key_type& finish) const
        {
            return key_range(basis, start, finish);
        }
    };

    TEST_FIXTURE(fixture, test_default_initialisation_equality)
    {
        iterator a, b;
        CHECK(a == b);
    }

    TEST_FIXTURE(fixture, test_default_initialisation_different_from_begin)
    {
        iterator a;
        auto basis_iter = iterate_basis();
        iterator b = basis_iter.begin();
        CHECK(a != b);
    }

    TEST_FIXTURE(fixture, test_begin_methods_yield_same_value)
    {
        auto basis_iter = iterate_basis();
        CHECK(basis_iter.begin() == basis_iter.begin());
    }

    TEST_FIXTURE(fixture, test_equal_after_increment)
    {
        auto basis_iter = iterate_basis();
        auto it1 = basis_iter.begin();
        auto it2 = basis_iter.begin();

        for (alg::DIMN i = 0; i < tensor_basis::start_of_degree(3); ++i) {
            CHECK(it1 == it2);
            ++it1;
            ++it2;
        }
    }

    TEST_FIXTURE(fixture, test_prefix_postfix_increment_equal)
    {
        auto basis_iter = iterate_basis();
        auto it1 = basis_iter.begin();
        auto it2 = basis_iter.begin();

        for (alg::DIMN i = 0; i < tensor_basis::start_of_degree(3); ++i) {
            CHECK(it1++ == it2);
            ++it2;
        }
    }

    TEST_FIXTURE(fixture, test_iterate_keys_deref_value)
    {
        auto key = basis.begin();

        for (auto r : iterate_basis()) {
            CHECK_EQUAL(key, r);
            key = basis.nextkey(key);
        }
    }

    TEST_FIXTURE(fixture, test_iteration_exhaustive)
    {
        alg::DIMN count = 0;
        for (auto _ UNUSED : iterate_basis()) {
            ++count;
        }
        CHECK_EQUAL(tensor_basis::start_of_degree(3), count);
    }

    TEST_FIXTURE(fixture, test_iteration_old_style)
    {
        auto key = basis.begin();
        auto basis_iter = iterate_basis();
        for (iterator it = basis_iter.begin(); it != basis_iter.end(); ++it) {
            CHECK_EQUAL(key, *it);
            key = basis.nextkey(key);
        }
    }

    TEST_FIXTURE(fixture, test_iteration_from_key)
    {
        key_type start_key{alg::LET(1), alg::LET(1)};// start of degree 2
        key_type k(start_key);
        auto basis_iter = iterate_basis(start_key);

        for (auto r : basis_iter) {
            CHECK_EQUAL(r, k);
            k = basis.nextkey(k);
        }
    }

    TEST_FIXTURE(fixture, test_iteration_between_keys)
    {
        key_type start_key{alg::LET(1)}, end_key{alg::LET{1}, alg::LET{1}};
        auto basis_iter = iterate_basis(start_key, end_key);

        auto first = basis_iter.begin();
        auto second = basis_iter.end();

        CHECK_EQUAL(start_key, *first);
        CHECK_EQUAL(end_key, *second);
        CHECK(++(++first) == second);

        key_type k{start_key};
        int count = 0;
        for (auto r : basis_iter) {
            ++count;
            CHECK_EQUAL(k, r);
            k = basis.nextkey(k);
        }
        CHECK_EQUAL(2, count);
    }
}
