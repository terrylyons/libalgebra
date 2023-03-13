
    TEST_FIXTURE(fixture, tensor_exp_zero) {
        TEST_DETAILS();
        TENSOR arg;
        TENSOR expected(S(1));

        CHECK_EQUAL(expected, exp(arg));
    }

    TEST_FIXTURE(fixture, tensor_exp_unit) {
        TEST_DETAILS();

        TENSOR arg(S(1));
        TENSOR expected(exp_to_depth(S(1), S(1)));

        CHECK_VEC_CLOSE(expected, exp(arg), expected_error);
    }

    TEST_FIXTURE(fixture, tensor_exponential_single_letter) {
        TEST_DETAILS();

        KEY key(KEY::LET(1UL));
        TENSOR arg(key);

        TENSOR expected(exp_to_depth(arg, S(1)));

        TENSOR result(exp(arg));

        CHECK_EQUAL(6, result.size());
        CHECK_VEC_CLOSE(expected, result, expected_error);
    }

    TEST_FIXTURE(fixture, tensor_exp_two_letters) {
        TEST_DETAILS();

        KEY key1(KEY::LET(1UL)), key2(KEY::LET(2UL));

        TENSOR arg;
        arg.add_scal_prod(key1, S(0.256147));
        arg.add_scal_prod(key2, S(0.852625));
        TENSOR result (exp(arg));

        TENSOR expected(exp_to_depth(arg, S(1)));

        // 63 = 1 + 2 + 4 + 8 + 16 + 32
        CHECK_EQUAL(63, result.size());
        CHECK_VEC_CLOSE(expected, result, expected_error);

    }

    TEST_FIXTURE(fixture, tensor_exp_two_letters_plus_unit) {
        TENSOR arg(S(1));
        KEY key1(KEY::LET(1UL)), key2(KEY::LET(2UL));
        arg.add_scal_prod(key1, S(0.256147));
        arg.add_scal_prod(key2, S(0.852625));
        TENSOR result (exp(arg));

        TENSOR expected(exp_to_depth(arg, S(1)));

        CHECK_EQUAL(63, result.size());
        CHECK_VEC_CLOSE(expected, result, expected_error);
    }


    TEST_FIXTURE(fixture, tensor_exp_degree_2_key) {
        TEST_DETAILS();

        KEY key(KEY::LET(1UL));
        key.push_back(KEY::LET(2UL));

        TENSOR arg;
        arg.add_scal_prod(key, S(0.534642));

        // exp should contain 1, (1,2), (1,2,1,2)
        TENSOR result(exp(arg));

        TENSOR expected(exp_to_depth(arg, S(1)));

        CHECK_EQUAL(3, result.size());
        CHECK_VEC_CLOSE(expected, result, expected_error);
    }

    TEST_FIXTURE(fixture, tensor_exp_degree_2_2_keys) {
        TEST_DETAILS();

        KEY key1(KEY::LET(1UL)), key2(KEY::LET(3UL));
        key1.push_back(KEY::LET(2UL));
        key2.push_back(KEY::LET(4UL));

        TENSOR arg;
        arg.add_scal_prod(key1, S(0.534642));
        arg.add_scal_prod(key2, S(0.163781));

        // exp should contain 1, (1,2), (3,4), (1,2,1,2), (1,2,3,4), (3,4,1,2), (3,4,3,4)
        TENSOR result(exp(arg));

        TENSOR expected(exp_to_depth(arg, S(1)));

        CHECK_EQUAL(7, result.size());
        CHECK_VEC_CLOSE(expected, result, expected_error);
    }


    // Tests for fused multiply exponentiate


    TEST_FIXTURE(fixture, fused_multiply_exp_zero_zero) {
        TENSOR t1, t2, zero;

        CHECK_EQUAL(zero, t1.fmexp(t2));
    }


    TEST_FIXTURE(fixture, fused_multiply_exp_one_zero) {
        KEY kunit;
        TENSOR t1(kunit), t2, expected(kunit);

        CHECK_EQUAL(expected, t1.fmexp(t2));
    }

    TEST_FIXTURE(fixture, fused_multiply_exp_equal_vs_old_single_key) {
        KEY kunit, k1(KEY::LET(1));
        TENSOR t1(kunit), t2(k1), expected(exp(t2));

        CHECK_EQUAL(expected, t1.fmexp(t2));
    }

    TEST_FIXTURE(fixture, fused_multiply_exp_equal_vs_old_two_keys) {
        KEY kunit, k1(KEY::LET(1)), k2(KEY::LET(2));
        TENSOR lhs(kunit), t1(k1), t2(k2), rhs(t1 + t2);

        CHECK_EQUAL(exp(rhs), lhs.fmexp(rhs));
    }

    TEST_FIXTURE(fixture, fused_multiply_exp_1_key_2_keys) {
        KEY klhs(KEY::LET(3)), k1(KEY::LET(1)), k2(KEY::LET(2));
        TENSOR lhs(klhs), t1(k1), t2(k2), rhs(t1 + t2);

        CHECK_EQUAL(lhs*exp(rhs), lhs.fmexp(rhs));
    }

    TEST_FIXTURE(fixture, fused_multiply_exp_1_key_2_keys_shared) {
        KEY klhs(KEY::LET(1)), k1(KEY::LET(1)), k2(KEY::LET(2));
        TENSOR lhs(klhs), t1(k1), t2(k2), rhs(t1 + t2);

        CHECK_EQUAL(lhs*exp(rhs), lhs.fmexp(rhs));
    }
