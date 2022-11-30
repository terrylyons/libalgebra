//
// Created by sam on 10/02/2021.
//

#include <UnitTest++.h>

#include <libalgebra/libalgebra.h>
#include <libalgebra/detail/integer_maths.h>

#include "../../common/time_and_details.h"


using alg::integer_maths::is_squarefree;
/*

SUITE(divisor_calculation) {

#define TEST_FIRST_DIVISOR(NUMBER, EXPECTED)                            \
    TEST(test_first_divisor_ ## NUMBER) {                               \
        TEST_DETAILS();                                                 \
        CHECK_EQUAL(EXPECTED, divisor_calc<NUMBER>::divisor);           \
        CHECK_EQUAL(NUMBER / EXPECTED, divisor_calc<NUMBER>::quotient); \
    }

    TEST_FIRST_DIVISOR(2, 2)
    TEST_FIRST_DIVISOR(3, 3)
    TEST_FIRST_DIVISOR(4, 2)
    TEST_FIRST_DIVISOR(5, 5)
    TEST_FIRST_DIVISOR(6, 2)
    TEST_FIRST_DIVISOR(7, 7)
    TEST_FIRST_DIVISOR(8, 2)
    TEST_FIRST_DIVISOR(9, 3)
    TEST_FIRST_DIVISOR(10, 2)
    TEST_FIRST_DIVISOR(11, 11)
    TEST_FIRST_DIVISOR(12, 2)
    TEST_FIRST_DIVISOR(13, 13)
    TEST_FIRST_DIVISOR(14, 2)
    TEST_FIRST_DIVISOR(15, 3)
    TEST_FIRST_DIVISOR(16, 2)
    TEST_FIRST_DIVISOR(17, 17)
    TEST_FIRST_DIVISOR(18, 2)
    TEST_FIRST_DIVISOR(19, 19)
    TEST_FIRST_DIVISOR(20, 2)
    TEST_FIRST_DIVISOR(21, 3)
    TEST_FIRST_DIVISOR(22, 2)
    TEST_FIRST_DIVISOR(23, 23)
    TEST_FIRST_DIVISOR(24, 2)
    TEST_FIRST_DIVISOR(25, 5)
    TEST_FIRST_DIVISOR(26, 2)


#undef TEST_FIRST_DIVISOR

}
*/

SUITE(square_free) {

#define TEST_SQUARE_FREE(NUMBER, EXPECTED)                                  \
    TEST(test_square_free_ ## NUMBER) {                                     \
        CHECK_EQUAL(EXPECTED, is_squarefree(NUMBER));                       \
    }
// there is a potential for error here unless square free returns bool does true == 3 ? No if true converts to an integer!
    TEST_SQUARE_FREE(2, true)
    TEST_SQUARE_FREE(3, true)
    TEST_SQUARE_FREE(4, false)
    TEST_SQUARE_FREE(5, true)
    TEST_SQUARE_FREE(6, true)
    TEST_SQUARE_FREE(7, true)
    TEST_SQUARE_FREE(8, false)
    TEST_SQUARE_FREE(9, false)
    TEST_SQUARE_FREE(10, true)
    TEST_SQUARE_FREE(11, true)
    TEST_SQUARE_FREE(12, false)
    TEST_SQUARE_FREE(13, true)
    TEST_SQUARE_FREE(14, true)
    TEST_SQUARE_FREE(15, true)
    TEST_SQUARE_FREE(16, false)
    TEST_SQUARE_FREE(17, true)
    TEST_SQUARE_FREE(18, false)
    TEST_SQUARE_FREE(19, true)
    TEST_SQUARE_FREE(20, false)
    TEST_SQUARE_FREE(21, true)
    TEST_SQUARE_FREE(22, true)
    TEST_SQUARE_FREE(23, true)
    TEST_SQUARE_FREE(24, false)
    TEST_SQUARE_FREE(25, false)
    TEST_SQUARE_FREE(26, true)
    TEST_SQUARE_FREE(27, false)
    TEST_SQUARE_FREE(28, false)

#undef TEST_SQUARE_FREE

}
