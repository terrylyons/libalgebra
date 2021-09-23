//
// Created by sam on 03/02/2021.
//

#if 0
#undef TEST_FIXTURE
#undef CHECK_EQUAL
#define CHECK_EQUAL CHECK_EQUAL
#define TEST_FIXTURE(F, N) void N()
#endif

#if defined(_VECTOR_TYPE)

    TEST_FIXTURE(Fixture, test_neutral_access_all_zero) {
        TEST_DETAILS();
        VECT neut;
        for (KEY i=0; i<BASIS::dimension; ++i)
            REQUIRE CHECK_EQUAL(S(0), neut[i]);
    }

    TEST_FIXTURE(Fixture, test_const_neutral_access_all_zero) {
        TEST_DETAILS();
        const VECT neut;

        for (KEY i=0; i<BASIS::dimension; ++i) {
            REQUIRE CHECK_EQUAL(S(0), neut[i]);
        }

    }

    TEST_FIXTURE(Fixture, test_unidim_set_element) {
        TEST_DETAILS();
        VECT vect(KEY(0));

        CHECK_EQUAL(S(1), vect[KEY(0)]);
        for (KEY i=1; i<BASIS::dimension; ++i)
            REQUIRE CHECK_EQUAL(S(0), vect[i]);
    }

    TEST_FIXTURE(Fixture, test_const_unidim_set_element) {
        TEST_DETAILS();
        const VECT vect(KEY(0));

        CHECK_EQUAL(S(1), vect[KEY(0)]);
        for (KEY i=1; i<BASIS::dimension; ++i)
            REQUIRE CHECK_EQUAL(S(0), vect[i]);
    }

#endif