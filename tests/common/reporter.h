//
// Created by sam on 26/03/2021.
//

#ifndef LIBALGEBRAUNITTESTS_REPORTER_H
#define LIBALGEBRAUNITTESTS_REPORTER_H
#include <iostream>

#include <UnitTest++.h>

class reporter : public UnitTest::TestReporter
{
public:
    ~reporter() override = default;

    void ReportTestStart(const UnitTest::TestDetails &test) override;
    void ReportTestFinish(const UnitTest::TestDetails &test, float secondsElapsed) override;
    void ReportFailure(const UnitTest::TestDetails &test, const char *failure) override;
    void ReportSummary(int totalTestCount, int failedTestCount, int failureCount, float secondsElapsed) override;

};


#endif //LIBALGEBRAUNITTESTS_REPORTER_H
