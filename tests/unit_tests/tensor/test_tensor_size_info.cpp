//
// Created by sam on 10/02/2021.
//

#include <UnitTest++.h>

#include <libalgebra/libalgebra.h>

#include "../../common/time_and_details.h"

SUITE(tensor_size_info) {

#define TEST_TENSOR_SIZE(WIDTH, DEPTH, EXPECTED)                \
    TEST(test_tensor_size_info_ ## WIDTH ## _ ## DEPTH) {       \
        TEST_DETAILS();                                         \
        alg::tensor_basis<WIDTH, DEPTH> basis;            \
        CHECK_EQUAL(EXPECTED, basis.start_of_degree(DEPTH));    \
    }

    TEST_TENSOR_SIZE(2, 0, 0)
    TEST_TENSOR_SIZE(2, 1, 1)
    TEST_TENSOR_SIZE(2, 2, 3)
    TEST_TENSOR_SIZE(2, 3, 7)
    TEST_TENSOR_SIZE(2, 4, 15)
    TEST_TENSOR_SIZE(2, 5, 31)
    TEST_TENSOR_SIZE(2, 6, 63)
    TEST_TENSOR_SIZE(2, 7, 127)
    TEST_TENSOR_SIZE(2, 8, 255)
    TEST_TENSOR_SIZE(2, 9, 511)
    TEST_TENSOR_SIZE(2, 10, 1023)
    TEST_TENSOR_SIZE(2, 11, 2047)
    TEST_TENSOR_SIZE(2, 12, 4095)
    TEST_TENSOR_SIZE(2, 13, 8191)
    TEST_TENSOR_SIZE(2, 14, 16383)
    TEST_TENSOR_SIZE(2, 15, 32767)
    TEST_TENSOR_SIZE(2, 16, 65535)
    TEST_TENSOR_SIZE(2, 17, 131071)
    TEST_TENSOR_SIZE(2, 18, 262143)
    TEST_TENSOR_SIZE(2, 19, 524287)
    TEST_TENSOR_SIZE(2, 20, 1048575)
    TEST_TENSOR_SIZE(2, 21, 2097151)
    TEST_TENSOR_SIZE(2, 22, 4194303)
    TEST_TENSOR_SIZE(2, 23, 8388607)
    TEST_TENSOR_SIZE(2, 24, 16777215)
    TEST_TENSOR_SIZE(2, 25, 33554431)
    TEST_TENSOR_SIZE(2, 26, 67108863)

    TEST_TENSOR_SIZE(3, 0, 0)
    TEST_TENSOR_SIZE(3, 1, 1)
    TEST_TENSOR_SIZE(3, 2, 4)
    TEST_TENSOR_SIZE(3, 3, 13)
    TEST_TENSOR_SIZE(3, 4, 40)
    TEST_TENSOR_SIZE(3, 5, 121)
    TEST_TENSOR_SIZE(3, 6, 364)
    TEST_TENSOR_SIZE(3, 7, 1093)
    TEST_TENSOR_SIZE(3, 8, 3280)
    TEST_TENSOR_SIZE(3, 9, 9841)
    TEST_TENSOR_SIZE(3, 10, 29524)
    TEST_TENSOR_SIZE(3, 11, 88573)
    TEST_TENSOR_SIZE(3, 12, 265720)
    TEST_TENSOR_SIZE(3, 13, 797161)
    TEST_TENSOR_SIZE(3, 14, 2391484)
    TEST_TENSOR_SIZE(3, 15, 7174453)
    TEST_TENSOR_SIZE(3, 16, 21523360)
    TEST_TENSOR_SIZE(3, 17, 64570081)

    TEST_TENSOR_SIZE(5, 0, 0)
    TEST_TENSOR_SIZE(5, 1, 1)
    TEST_TENSOR_SIZE(5, 2, 6)
    TEST_TENSOR_SIZE(5, 3, 31)
    TEST_TENSOR_SIZE(5, 4, 156)
    TEST_TENSOR_SIZE(5, 5, 781)
    TEST_TENSOR_SIZE(5, 6, 3906)
    TEST_TENSOR_SIZE(5, 7, 19531)
    TEST_TENSOR_SIZE(5, 8, 97656)
    TEST_TENSOR_SIZE(5, 9, 488281)
    TEST_TENSOR_SIZE(5, 10, 2441406)


#undef TEST_TENSOR_SIZE

}
