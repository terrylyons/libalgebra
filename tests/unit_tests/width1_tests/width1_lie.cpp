//
// Created by user on 20/04/23.
//

#include "width1.h"

SUITE(Width1)
{
    TEST_FIXTURE(Width1Tests, test_lie_key_creation)
    {
        typename LIE::KEY lkey(1);

        LIE lie(lkey);

        CHECK_EQUAL(1.0, lie[lkey]);
    }

    TEST_FIXTURE(Width1Tests, test_hall_set_size) {

        CHECK_EQUAL(1, LIE::basis.size());

    }

}
