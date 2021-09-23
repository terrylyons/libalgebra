//
// Created by sam on 12/02/2021.
//

#ifdef ALGEBRA_TESTS_VECT_TYPE



TEST_FIXTURE(Fixture, test_zero_zero_mul) {
    TEST_DETAILS();
    ALG lhs, rhs, expected;

    CHECK_EQUAL(expected, lhs*rhs);

}

TEST_FIXTURE(Fixture, test_unidim_same_keys) {
    TEST_DETAILS();
    S s1 = rand_scalar(), s2 = rand_scalar();
    KEY k1 = 1, k2 = 1;
    ALG lhs(k1, s1), rhs(k2, s2), expected(k1, s1*s2);

    CHECK_EQUAL(expected, lhs*rhs);
}

TEST_FIXTURE(Fixture, test_unidim_different_keys) {
    TEST_DETAILS();
    S s1 = rand_scalar(), s2 = rand_scalar();
    KEY k1 = 0, k2 = 1;
    ALG lhs(k1, s1), rhs(k2, s2), expected;

    CHECK_EQUAL(expected, lhs*rhs);
}

TEST_FIXTURE(Fixture, test_full_random_alg) {
    TEST_DETAILS();
    ALG lhs = rand_vec(), rhs = rand_vec(), expected;
    for (KEY i=0; i < BASIS::dimension; ++i) {
        expected[i] = lhs[i] * rhs[i];
    }

    CHECK_EQUAL(expected, lhs*rhs);
}

TEST_FIXTURE(Fixture, test_inplace_equals_external) {
    TEST_DETAILS();

    ALG lhs = rand_vec(), rhs = rand_vec(), expected(lhs);

    expected *= rhs;
    CHECK_EQUAL(expected, lhs * rhs);
}

TEST_FIXTURE(Fixture, test_fused_add_mul) {
    TEST_DETAILS();

    ALG a = rand_vec(), b = rand_vec(), c = rand_vec();
    ALG expected = a + (b * c);

    CHECK_EQUAL(expected, a.add_mul(b, c));
}

TEST_FIXTURE(Fixture, test_fused_sub_mul) {
    TEST_DETAILS();

    ALG a = rand_vec(), b = rand_vec(), c = rand_vec();
    ALG expected = a - (b * c);

    CHECK_EQUAL(expected, a.sub_mul(b, c));
}

TEST_FIXTURE(Fixture, test_mul_scal_prod) {
    TEST_DETAILS();
    ALG lhs = rand_vec(), rhs = rand_vec();
    S scal = rand_scalar(S(1), S(5)); // non-zero

    ALG expected = lhs * (rhs * scal);

    CHECK_EQUAL(expected, lhs.mul_scal_prod(rhs, scal));
}

TEST_FIXTURE(Fixture, test_mul_scal_div) {
    TEST_DETAILS();
    ALG lhs = rand_vec(), rhs = rand_vec();
    Q scal = rand_scalar(S(1), S(5)); // non-zero

    ALG expected = lhs * (rhs / scal);
    CHECK_EQUAL(expected, lhs.mul_scal_div(rhs, scal));
}

TEST_FIXTURE(Fixture, test_commutator) {
    TEST_DETAILS();
    ALG a = rand_vec(), b = rand_vec();

    ALG expected = a * b - b * a;

    CHECK_EQUAL(expected, commutator(a, b));

}

#endif