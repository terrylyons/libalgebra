//
// Created by sam on 22/02/2021.
//

#if defined(_VECTOR_TYPE)


TEST_FIXTURE(Fixture, test_L1_norm_neutral_0) {
    TEST_DETAILS();

    VECT neut;

    CHECK_EQUAL(S(0), neut.NormL1());
}


TEST_FIXTURE(Fixture, test_L1_norm_unidim_positive_scalar) {
    TEST_DETAILS();

    S rval (rand_scalar(1, 5));
    VECT vect(KEY(1), rval);

    CHECK_EQUAL(rval, vect.NormL1());
}

TEST_FIXTURE(Fixture, test_L1_norm_unidim_negative_scalar) {
    TEST_DETAILS();

    S rval(rand_scalar(1, 5));
    VECT vect(KEY(1), -rval);

    CHECK_EQUAL(rval, vect.NormL1());
}


TEST_FIXTURE(Fixture, test_L1_norm_unidim_set_0) {
    TEST_DETAILS();

    VECT vect(KEY(1), S(0));

    CHECK_EQUAL(S(0), vect.NormL1());
}

TEST_FIXTURE(Fixture, test_L1_norm_random_vect) {
    TEST_DETAILS();
    VECT vect(rand_vec());

    S ans(0);
    for (KEY i(0); i<BASIS::dimension; ++i) {
        ans += abs(vect[i]);
    }

    CHECK_EQUAL(ans, vect.NormL1());
}



TEST_FIXTURE(Fixture, test_Linf_norm_neutral_0) {
    TEST_DETAILS();

    VECT neut;

    CHECK_EQUAL(S(0), neut.NormLInf());
}


TEST_FIXTURE(Fixture, test_Linf_norm_unidim_positive_scalar) {
    TEST_DETAILS();

    S rval (rand_scalar(1, 5));
    VECT vect(KEY(1), rval);

    CHECK_EQUAL(rval, vect.NormLInf());
}

TEST_FIXTURE(Fixture, test_Linf_norm_unidim_negative_scalar) {
    TEST_DETAILS();

    S rval(rand_scalar(1, 5));
    VECT vect(KEY(1), -rval);

    CHECK_EQUAL(rval, vect.NormLInf());
}


TEST_FIXTURE(Fixture, test_Linf_norm_unidim_set_0) {
    TEST_DETAILS();

    VECT vect(KEY(1), S(0));

    CHECK_EQUAL(S(0), vect.NormLInf());
}

TEST_FIXTURE(Fixture, test_Linf_norm_random_vect) {
    TEST_DETAILS();
    VECT vect(rand_vec());

    S ans(0);
    for (KEY i(0); i<BASIS::dimension; ++i) {
        ans = std::max(ans, abs(vect[i]));
    }

    CHECK_EQUAL(ans, vect.NormLInf());
}





#endif //LIBALGEBRAUNITTESTS_VECTOR_NORM_SUITE_H
