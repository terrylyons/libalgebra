//
// Created by sam on 02/02/2021.
//

#include <UnitTest++.h>
#include "reporter.h"

#include <vector>
#include <cstring>

struct FilterList : std::vector<std::pair<char*, char*>> {

    ~FilterList() {
        for (auto&& item : *this) {
            delete[] item.first;
            delete[] item.second;
        }
    }

};

static void add_filter(FilterList& filters, char* text) {
    char* split = strchr(text, ':');
    char *suite = nullptr, *name = nullptr;

    if (split != nullptr) {
        ++split;// Fine since split starts at ':'
        auto sizeof_suite = static_cast<size_t>(split - text);
        auto sizeof_name = strlen(split);

        suite = new char[sizeof_suite + 1];
        name = new char[sizeof_name + 1];

        strcpy(name, split);
        if (sizeof_suite > 1) {
            memcpy(suite, text, sizeof_suite - 1);
        }
        suite[sizeof_suite] = '\0';
    }
    else {
        auto len = strlen(text);
        name = new char[len + 1];
        strcpy(name, text);
    }

    if (name != nullptr || suite != nullptr) {
        filters.emplace_back(suite, name);
    }
}

int main(int argc, char** argv)
{
    FilterList filters;

    for (int i=1; i<argc; ++i) {
        char* arg = argv[i];
        if (strncmp(arg, "--test-filter", 13) == 0) {
            if (arg[13] == '\0') {
                if ((i + 1) >= argc || strncmp(argv[i + 1], "-", 1) == 0) {
                    std::cout << "no argument provided to --filter\n";
                    return 1;
                }
                add_filter(filters, argv[i + 1]);
            } else if (arg[13] == '=') {
                char* real_arg = arg + 14;// "--test-filter=" is 13 characters, so 14 is safe
                if (*real_arg == '\0') {
                    std::cout << "Expected non-empty string after =\n";
                    return 1;
                }
                add_filter(filters, real_arg);

            } else {
                std::cout << "Unrecognised filter format\n";
                return 1;
            }


        } else if (strcmp(arg, "--list-tests") == 0) {
            auto test_list = UnitTest::Test::GetTestList();

            auto* test = test_list.GetHead();

            while (test != nullptr) {
                std::cout << test->m_details.suiteName
                          << ':'
                          << test->m_details.testName
                          << '\n';

                test = test->m_nextTest;
            }
            return 0;
        }
    }


    reporter report;
    UnitTest::TestRunner runner(report);

    auto predicate = [&filters](UnitTest::Test* test) {
        for (const auto& filter : filters) {
            if (strcmp(filter.first, test->m_details.suiteName) == 0 && strcmp(filter.second, test->m_details.testName) == 0) {
                return true;
            }
        }
        return false;
    };

    return runner.RunTestsIf(UnitTest::Test::GetTestList(), NULL, predicate, 0);

}
