//
// Created by sam on 17/11/2021.
//

#ifndef LIBALGEBRA_MULTI_TEST_H
#define LIBALGEBRA_MULTI_TEST_H

#include <UnitTest++.h>
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

struct metadata {
    int lineno;
    const char* currentName;
    const char* currentSuite;
    const char* currentFile;
};

using test_pair = std::pair<metadata, std::unique_ptr<UnitTest::Test>>;

class generative_suite
{

    std::vector<test_pair>& tests;

    metadata current;

    template<typename Fn>
    struct test_type : public UnitTest::Test {
        explicit test_type(const metadata& m, Fn&& fn)
            : UnitTest::Test(m.currentName, m.currentSuite, m.currentFile, m.lineno),
              m_impl(std::forward<Fn>(fn))
        {}

        void RunImpl() const override
        {
            m_impl();
        }

    private:
        Fn m_impl;
    };

public:
    explicit generative_suite(std::vector<test_pair>& t, const char* suite) : tests(t), current{0, nullptr, suite, nullptr}
    {}

    void set_metadata(const char* name, const char* filename, int lineno)
    {
        current.currentName = name;
        current.currentFile = filename;
        current.lineno = lineno;
    }

    template<typename Fn>
    void operator+(Fn&& fn)
    {
        tests.push_back(
                std::make_pair(
                        current,
                        std::unique_ptr<UnitTest::Test>(new test_type<Fn>(current, std::forward<Fn>(fn)))));
    }

    ~generative_suite()
    {
        auto& list = UnitTest::Test::GetTestList();
        for (auto& t : tests) {

            list.Add(t.second.get());
        }
    }
};

#define NEW_AUTO_SUITE(NAME, ...)                                          \
    template<__VA_ARGS__>                                                  \
    void create_tests_##NAME(la_testing::generative_suite& s);             \
                                                                           \
    template<typename... Args>                                             \
    std::vector<la_testing::test_pair> make_suite_##NAME(const char* name) \
    {                                                                      \
        std::vector<la_testing::test_pair> tests;                          \
                                                                           \
        la_testing::generative_suite suite(tests, name);                   \
        create_tests_##NAME<Args...>(suite);                               \
                                                                           \
        return tests;                                                      \
    }                                                                      \
                                                                           \
    template<__VA_ARGS__>                                                  \
    void create_tests_##NAME(la_testing::generative_suite& s)

#define ADD_TEST(NAME)                         \
    s.set_metadata(#NAME, __FILE__, __LINE__); \
    s + []()

#define MAKE_SUITE_FOR(INSTANCE_NAME, SUITE_NAME, ...) \
    static auto INSTANCE_NAME = make_suite_##SUITE_NAME<__VA_ARGS__>(#INSTANCE_NAME);

}// namespace la_testing

#endif//LIBALGEBRA_MULTI_TEST_H
