    TEST_FIXTURE(fixture, tensor_log_tensor_unit) {
        TEST_DETAILS();

        KEY kunit;
        TENSOR tunit(kunit);

        TENSOR expected;

        CHECK_EQUAL(expected, log(tunit));
    }

    TEST_FIXTURE(fixture, tensor_log_zero_tensor) {
        TEST_DETAILS();

        // Implicitly, this is not the log of 0, but of 1 + 0;
        TENSOR arg, expected;

        CHECK_EQUAL(expected, log(arg));
    }

    TEST_FIXTURE(fixture, tensor_log_single_letter) {
        TEST_DETAILS();

        TENSOR arg(1UL, S(1));

        TENSOR result(log(arg));

        TENSOR expected(log_to_depth(arg));

        CHECK_EQUAL(5, result.size());
        CHECK_VEC_CLOSE(expected, result, expected_error);
    }

    TEST_FIXTURE(fixture, tensor_log_2_keys) {
        TEST_DETAILS();

        KEY key1(KEY::LET(1UL)), key2(KEY::LET(2UL));
        TENSOR arg;
        arg.add_scal_prod(key1, S(0.647183));
        arg.add_scal_prod(key2, S(0.368154));

        TENSOR result(log(arg));
        TENSOR expected(log_to_depth(arg));

        // 62 = 2 + 4 + 8 + 16 + 32
        CHECK_EQUAL(62, result.size());
        CHECK_VEC_CLOSE(expected, result, expected_error);
    }

    TEST_FIXTURE(fixture, tensor_log_2_keys_and_unit) {
        TEST_DETAILS();

        KEY key1(KEY::LET(1UL)), key2(KEY::LET(2UL)), kunit;
        TENSOR arg;
        arg.add_scal_prod(kunit, S(1));
        arg.add_scal_prod(key1, S(0.647183));
        arg.add_scal_prod(key2, S(0.368154));

        TENSOR result(log(arg));
        TENSOR expected(log_to_depth(arg));

        // 62 = 2 + 4 + 8 + 16 + 32
        CHECK_EQUAL(62, result.size());
        CHECK_VEC_CLOSE(expected, result, expected_error);
    }

    TEST_FIXTURE(fixture, tensor_log_degree_2_key) {
        TEST_DETAILS();

        KEY key;
        key.push_back(KEY::LET(1UL));
        key.push_back(KEY::LET(2UL));

        TENSOR arg;
        arg.add_scal_prod(key, S(0.212346));

        TENSOR result(log(arg));
        TENSOR expected(log_to_depth(arg));

        // Log should contain (1,2) (1,2,1,2)
        CHECK_EQUAL(2, result.size());
        CHECK_VEC_CLOSE(expected, result, expected_error);
    }

    TEST_FIXTURE(fixture, tensor_log_degree_2_key_and_unit) {
        TEST_DETAILS();

        KEY key, kunit;
        key.push_back(KEY::LET(1UL));
        key.push_back(KEY::LET(2UL));

        TENSOR arg;
        arg.add_scal_prod(kunit, S(1));
        arg.add_scal_prod(key, S(0.212346));

        TENSOR result(log(arg));
        TENSOR expected(log_to_depth(arg));

        // Log should contain (1,2) (1,2,1,2)
        CHECK_EQUAL(2, result.size());
        CHECK_VEC_CLOSE(expected, result, expected_error);
    }

    TEST_FIXTURE(fixture, tensor_log_degree_2_2_keys) {
        TEST_DETAILS();

        KEY key1, key2;
        key1.push_back(KEY::LET(1UL));
        key1.push_back(KEY::LET(2UL));
        key2.push_back(KEY::LET(3UL));
        key2.push_back(KEY::LET(4UL));

        TENSOR arg;
        arg.add_scal_prod(key1, S(0.898525));
        arg.add_scal_prod(key2, S(0.613454));

        TENSOR result(log(arg));
        TENSOR expected(log_to_depth(arg));

        // Result should contain:
        // (1,2), (3,4), (1,2,1,2), (1,2,3,4), (3,4,1,2), (3,4,3,4)
        CHECK_EQUAL(6, result.size());
        CHECK_VEC_CLOSE(expected, result, expected_error);
    }

    TEST_FIXTURE(fixture, tensor_log_degree_2_2_keys_and_unit) {
        TEST_DETAILS();

        KEY key1, key2, kunit;
        key1.push_back(KEY::LET(1UL));
        key1.push_back(KEY::LET(2UL));
        key2.push_back(KEY::LET(3UL));
        key2.push_back(KEY::LET(4UL));

        TENSOR arg;
        arg.add_scal_prod(kunit, S(1));
        arg.add_scal_prod(key1, S(0.898525));
        arg.add_scal_prod(key2, S(0.613454));

        TENSOR result(log(arg));
        TENSOR expected(log_to_depth(arg));

        // Result should contain:
        // (1,2), (3,4), (1,2,1,2), (1,2,3,4), (3,4,1,2), (3,4,3,4)
        CHECK_EQUAL(6, result.size());
        CHECK_VEC_CLOSE(expected, result, expected_error);
    }
