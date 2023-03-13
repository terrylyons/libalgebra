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

    TEST_FIXTURE(Fixture, test_neutral_elements_equal) {
        TEST_DETAILS();
        VECT lhs, rhs;

        CHECK_EQUAL(lhs, rhs);
    }

    TEST_FIXTURE(Fixture, test_neutral_element_not_equal_to_non_neutral) {
        TEST_DETAILS();
        VECT lhs, rhs = rand_vec();
        rhs[0] = S(1);

        CHECK(!(lhs == rhs));
    }

    TEST_FIXTURE(Fixture, test_copy_constructed_elements_test_equal) {
        TEST_DETAILS();
        VECT lhs = rand_vec(), rhs(lhs);

        CHECK_EQUAL(lhs, rhs);
    }

    TEST_FIXTURE(Fixture, test_unidim_different_keys_not_equal) {
        TEST_DETAILS();
        VECT lhs(KEY(0)), rhs(KEY(1));

        CHECK(!(lhs == rhs));
    }

    TEST_FIXTURE(Fixture, test_unidim_same_key_different_values_not_equal) {
        TEST_DETAILS();
        VECT lhs(KEY(0), S(1)), rhs(KEY(0), S(1.5));

        CHECK(!(lhs == rhs));
    }


#endif