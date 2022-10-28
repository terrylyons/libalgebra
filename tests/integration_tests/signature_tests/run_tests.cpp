//
// Created by sam on 02/02/2021.
//

#include <UnitTest++.h>
#include "../../common/reporter.h"

int main()
{
    reporter report;
    UnitTest::TestRunner runner(report);

    return runner.RunTestsIf(UnitTest::Test::GetTestList(), NULL, UnitTest::True(), 0);
}
