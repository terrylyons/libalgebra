//
// Created by sam on 26/03/2021.
//

#include "reporter.h"
#include <cmath>


void reporter::ReportTestStart(const UnitTest::TestDetails &test)
{
    std::cout << "\nFile : " << test.filename
              << "\nLine : " << test.lineNumber
              << "\nSuite : " << test.suiteName
              << "\nTest : " << test.testName
              << std::endl;
}


void reporter::ReportFailure(const UnitTest::TestDetails &test, const char *failure)
{
    using namespace std;
#if defined(__APPLE__) || defined(__GNUG__)
    char const* const errorFormat = "%s:%d:%d: error: Failure in %s: %s\n";
    fprintf(stderr, errorFormat, test.suiteName, test.lineNumber, 1, test.testName, failure);
#else
    char const* const errorFormat = "%s(%d): error: Failure in %s: %s\n";
      fprintf(stderr, errorFormat, test.filename, test.lineNumber, test.testName, failure);
#endif
}

void reporter::ReportTestFinish(const UnitTest::TestDetails &test, float secondsElapsed)
{

    std::cout << "Elapsed time: " << secondsElapsed << " seconds" << std::endl;
}

void reporter::ReportSummary(int totalTestCount, int failedTestCount, int failureCount, float secondsElapsed)
{
    using namespace std;

    if (failureCount > 0)
        printf("FAILURE: %d out of %d tests failed (%d failures).\n", failedTestCount, totalTestCount, failureCount);
    else
        printf("Success: %d tests passed.\n", totalTestCount);

    printf("Test time: %.2f seconds.\n", secondsElapsed);
}
