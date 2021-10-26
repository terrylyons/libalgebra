TEST_FIXTURE(fixture, tensor_inverse_tunit) {
TEST_DETAILS();

KEY kunit;
TENSOR tunit(kunit);

TENSOR result(inverse(tunit));

CHECK_EQUAL(tunit, result);
}

TEST_FIXTURE(fixture, test_inverse_single_letter) {
TEST_DETAILS();

KEY kunit, key(KEY::LET(1UL));

TENSOR arg(kunit), tunit(kunit);
arg.add_scal_prod(key, S(0.264735));


CHECK_VEC_CLOSE(tunit, inverse(arg) * arg, expected_error);
}

TEST_FIXTURE(fixture, test_inverse_2_letters) {
TEST_DETAILS();

KEY kunit, key1(KEY::LET(1UL)), key2(KEY::LET(2UL));

TENSOR arg(kunit), tunit(kunit);
arg.add_scal_prod(key1, S(0.264735));
arg.add_scal_prod(key2, S(0.783473));

CHECK_VEC_CLOSE(tunit, inverse(arg) * arg, expected_error);
}

TEST_FIXTURE(fixture, test_inverse_degree_1_and_2_keys) {
TEST_DETAILS();

KEY kunit, key1(KEY::LET(1UL)), key2(KEY::LET(2UL));
TENSOR arg(kunit), tunit(kunit);
arg.add_scal_prod(key1, S(0.846539));
arg.add_scal_prod(key2, S(0.724562));

CHECK_VEC_CLOSE(tunit, inverse(arg) * arg, expected_error);
}

TEST_FIXTURE(fixture, test_inverse_involution_letter) {
TEST_DETAILS();

KEY kunit, key(KEY::LET(1UL));
TENSOR arg(kunit);
arg.add_scal_prod(key, S(0.362374));

CHECK_VEC_CLOSE(arg, inverse(inverse(arg)), expected_error);
}

TEST_FIXTURE(fixture, test_inverse_involution_deg_2_key) {
TEST_DETAILS();

KEY kunit, key(KEY::LET(1UL));
key.push_back(KEY::LET(2UL));
TENSOR arg(kunit);
arg.add_scal_prod(key, S(0.362374));

CHECK_VEC_CLOSE(arg, inverse(inverse(arg)), expected_error);
}

TEST_FIXTURE(fixture, test_inverse_involution_deg_1_and_2_keys) {
TEST_DETAILS();

KEY kunit, key1(KEY::LET(1UL)), key2(KEY::LET(2UL));
key2.push_back(KEY::LET(3UL));
TENSOR arg(kunit);
arg.add_scal_prod(key1, S(0.362374));
arg.add_scal_prod(key2, S(0.613471));

CHECK_VEC_CLOSE(arg, inverse(inverse(arg)), expected_error);
}
