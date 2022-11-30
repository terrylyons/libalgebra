//
// Created by sam on 17/11/2021.
//

#ifndef LIBALGEBRA_BASIS_ITERATION_SUITE_H
#define LIBALGEBRA_BASIS_ITERATION_SUITE_H

#include <UnitTest++.h>

#include <multi_test.h>
#include <string>

template<typename Basis>
la_testing::multi_suite create_basis_test_suite(const std::string& basis_name)
{
    la_testing::multi_suite suite;

    using key_type = typename Basis::KEY;

    suite.emplace_back(la_testing::make_multi_test(
            "test_" + basis_name + "_basic_iteration",
            basis_name + "_iteration_tests", __FILE__, __LINE__,
            []() {
                Basis basis;

                key_type k = basis.begin();
                for (auto r : basis.iterate_keys()) {
                    CHECK_EQUAL(k, r);
                    k = basis.nextkey(k);
                }
            }));

    suite.emplace_back(la_testing::make_multi_test(
            "test_" + basis_name + "_between_iteration",
            basis_name + "_iteration_tests", __FILE__, __LINE__,
            []() {
                Basis basis;
                key_type start(basis.index_to_key(1U)), finish(basis.index_to_key(3U));

                key_type k(start);
                for (auto r : basis.iterate_keys(start, finish)) {
                    CHECK_EQUAL(k, r);
                    k = basis.nextkey(k);
                }
                CHECK_EQUAL(k, finish);
            }));

    auto& list = UnitTest::Test::GetTestList();
    for (auto& tc : suite) {
        list.Add(tc.get());
    }

    return suite;
}

NEW_AUTO_SUITE(basis_iteration_tests, typename Basis)
{
    using key_type = typename Basis::KEY;

    ADD_TEST(test_basic_iteration) {
        Basis basis;

        key_type k = basis.begin();
        for (auto r : basis.iterate_keys()) {
            CHECK_EQUAL(k, r);
            k = basis.nextkey(k);
        }
    };


    ADD_TEST(test_between_iteration) {
        Basis basis;
        key_type start(basis.index_to_key(1U)), finish(basis.index_to_key(3U));

        key_type k(start);
        for (auto r : basis.iterate_keys(start, finish)) {
            CHECK_EQUAL(k, r);
            k = basis.nextkey(k);
        }
        CHECK_EQUAL(k, finish);
    };


}

#endif//LIBALGEBRA_BASIS_ITERATION_SUITE_H
