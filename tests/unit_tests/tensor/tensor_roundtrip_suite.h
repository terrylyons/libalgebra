TEST_FIXTURE(fixture, tensor_log_exp_roundtrip_letter) {
TEST_DETAILS();

KEY key(KEY::LET(1UL));

TENSOR arg(key);

TENSOR result(log(exp(arg)));

CHECK_VEC_CLOSE(arg, result, expected_error);
}

TEST_FIXTURE(fixture, tensor_log_exp_roundtrip_2_letter) {
TEST_DETAILS();

KEY key1(KEY::LET(1UL)), key2(KEY::LET(2UL));

TENSOR arg;
arg.add_scal_prod(key1, S(0.627431));
arg.add_scal_prod(key2, S(0.725713));

TENSOR result(log(exp(arg)));

CHECK_VEC_CLOSE(arg, result, expected_error);
}

TEST_FIXTURE(fixture, tensor_log_exp_roundtrip_degree_2_key) {
TEST_DETAILS();

KEY key(KEY::LET(1UL));
key.push_back(KEY::LET(2UL));

TENSOR arg;
arg.add_scal_prod(key, S(0.627431));

TENSOR result(log(exp(arg)));

CHECK_VEC_CLOSE(arg, result, expected_error);
}

TEST_FIXTURE(fixture, tensor_log_exp_roundtrip_deg_1_and_2_keys) {
TEST_DETAILS();

KEY key1(KEY::LET(1UL)), key2(KEY::LET(2UL));
key2.push_back(KEY::LET(3UL));

TENSOR arg;
arg.add_scal_prod(key1, S(0.724582));
arg.add_scal_prod(key2, S(0.623783));

TENSOR result(log(exp(arg)));

CHECK_VEC_CLOSE(arg, result, expected_error);
}

TEST_FIXTURE(fixture, tensor_log_exp_roundtrip_degree_2_2_keys) {
TEST_DETAILS();

KEY key1(KEY::LET(1UL)), key2(KEY::LET(3UL));
key1.push_back(KEY::LET(2UL));
key2.push_back(KEY::LET(4UL));

TENSOR arg;
arg.add_scal_prod(key1, S(0.627431));
arg.add_scal_prod(key2, S(0.735234));

TENSOR result(log(exp(arg)));

CHECK_VEC_CLOSE(arg, result, expected_error);
}


TEST_FIXTURE(fixture, tensor_exp_log_roundtrip_letter) {
TEST_DETAILS();

KEY key(KEY::LET(1UL)), kunit;

TENSOR arg(key), tunit(kunit);

TENSOR result(exp(log(arg)));

CHECK_VEC_CLOSE(arg + tunit, result, expected_error);
}

TEST_FIXTURE(fixture, tensor_exp_log_roundtrip_2_letter) {
TEST_DETAILS();

KEY key1(KEY::LET(1UL)), key2(KEY::LET(2UL)), kunit;

TENSOR arg, tunit(kunit);
arg.add_scal_prod(key1, S(0.627431));
arg.add_scal_prod(key2, S(0.725713));

TENSOR result(exp(log(arg)));

CHECK_VEC_CLOSE(arg + tunit, result, expected_error);
}

TEST_FIXTURE(fixture, tensor_exp_log_roundtrip_degree_2_key) {
TEST_DETAILS();

KEY key(KEY::LET(1UL)), kunit;
key.push_back(KEY::LET(2UL));

TENSOR arg, tunit(kunit);
arg.add_scal_prod(key, S(0.627431));

TENSOR result(exp(log(arg)));

CHECK_VEC_CLOSE(arg + tunit, result, expected_error);
}

TEST_FIXTURE(fixture, tensor_exp_log_roundtrip_deg_1_and_2_keys) {
TEST_DETAILS();

KEY key1(KEY::LET(1UL)), key2(KEY::LET(2UL)), kunit;
key2.push_back(KEY::LET(3UL));

TENSOR arg, tunit(kunit);
arg.add_scal_prod(key1, S(0.724582));
arg.add_scal_prod(key2, S(0.623783));

TENSOR result(exp(log(arg)));

CHECK_VEC_CLOSE(arg + tunit, result, expected_error);
}

TEST_FIXTURE(fixture, tensor_exp_log_roundtrip_degree_2_2_keys) {
TEST_DETAILS();

KEY key1(KEY::LET(1UL)), key2(KEY::LET(3UL)), kunit;
key1.push_back(KEY::LET(2UL));
key2.push_back(KEY::LET(4UL));


TENSOR arg, tunit(kunit);
arg.add_scal_prod(key1, S(0.627431));
arg.add_scal_prod(key2, S(0.735234));

TENSOR result(exp(log(arg)));

CHECK_VEC_CLOSE(arg + tunit, result, expected_error);
}