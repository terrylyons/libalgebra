//
// Created by sam on 17/11/2021.
//

#ifndef LIBALGEBRA_MULTI_TEST_H
#define LIBALGEBRA_MULTI_TEST_H

#include <UnitTest++/UnitTest++.h>
#include <memory>
#include <string>
#include <vector>

namespace la_testing {

template<typename TestImpl>
class multi_test : public UnitTest::Test
{
    TestImpl m_impl;

    std::string m_name;
    std::string m_suite;
    std::string m_file;

public:
    explicit multi_test(std::string&& testName, std::string&& suiteName, std::string&& filename, int lineNumber, TestImpl impl)
        : UnitTest::Test(testName.c_str(), suiteName.c_str(), filename.c_str(), lineNumber), m_impl(impl), m_name(std::move(testName)), m_suite(std::move(suiteName)),
          m_file(std::move(filename))
    {}

    void RunImpl() const override
    {
        m_impl();
    }
};

template<typename Impl>
std::unique_ptr<UnitTest::Test> make_multi_test(std::string&& testName, std::string&& suiteName, std::string&& filename, int lineno, Impl impl)
{
    return std::unique_ptr<UnitTest::Test>(new multi_test<Impl>(std::move(testName), std::move(suiteName), std::move(filename), lineno, impl));
}

using multi_suite = std::vector<std::unique_ptr<UnitTest::Test>>;
}// namespace la_testing

#endif//LIBALGEBRA_MULTI_TEST_H
