// LibAlgebraUnitTests.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

  // test.cpp
#include <UnitTest++.h>
// file:///C:/Program%20Files/dev/unittest-cpp-1.4-v90-v140/unittest-cpp-1.4/UnitTest++/docs/UnitTest++.html
//TEST(FailSpectacularly)
//{
//	CHECK(false);
//}
#include "../common/reporter.h"

int main()
{
    reporter report;
    UnitTest::TestRunner runner(report);

    return runner.RunTestsIf(UnitTest::Test::GetTestList(), NULL, UnitTest::True(), 0);
}
