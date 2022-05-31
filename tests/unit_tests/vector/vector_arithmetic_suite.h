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

#define TEST_BINARY_NEUTRAL(NAME, OP)                            \
    TEST_FIXTURE(Fixture, NAME ## _neutral_neutral) {            \
        TEST_DETAILS();                                          \
        VECT lhs, rhs, expected;                                 \
                                                                 \
        CHECK_EQUAL(expected, lhs OP rhs);                       \
    }                                                            \
                                                                 \
    TEST_FIXTURE(Fixture, NAME ## _random_neutral) {             \
        TEST_DETAILS();                                          \
        VECT neut, other = rand_vec(), expected(other);          \
                                                                 \
        CHECK_EQUAL(expected, other OP neut);                    \
    }                                                            \
                                                                 \
    TEST_FIXTURE(Fixture, NAME ## _netural_random) {             \
        TEST_DETAILS();                                          \
        VECT neut, other = rand_vec(), expected;                 \
        for (KEY i=0; i<BASIS::dimension; ++i)                   \
            expected[i] = S(0) OP other[i];                      \
                                                                 \
        CHECK_EQUAL(expected, neut OP other);                    \
    }

TEST_BINARY_NEUTRAL(test_addition, +)

TEST_BINARY_NEUTRAL(test_subtraction, -)


#undef TEST_BINARY_NEUTRAL

#define TEST_BINARY_OPERATOR(NAME, OP)                           \
    TEST_FIXTURE(Fixture, NAME ## _random_random) {              \
        TEST_DETAILS();                                          \
        VECT lhs = rand_vec(), rhs = rand_vec(), expected;       \
        for (KEY i=0; i < BASIS::dimension; ++i)                 \
            expected[i] = lhs[i] OP rhs[i];                      \
                                                                 \
        CHECK_EQUAL(expected, lhs OP rhs);                       \
    }

TEST_BINARY_OPERATOR(test_addition, +)

TEST_BINARY_OPERATOR(test_subtraction, -)

#undef TEST_BINARY_OPERATOR

#define TEST_INPLACE_BINARY_NEUTRAL(NAME, OP)                    \
    TEST_FIXTURE(Fixture, NAME ## _neutral_neutral) {            \
        TEST_DETAILS();                                          \
        VECT lhs, rhs, expected;                                 \
        lhs OP rhs;                                              \
        CHECK_EQUAL(expected, lhs);                              \
    }                                                            \
                                                                 \
    TEST_FIXTURE(Fixture, NAME ## _random_neutral) {             \
        TEST_DETAILS();                                          \
        VECT neut, other = rand_vec(), expected(other);          \
        other OP neut;                                           \
                                                                 \
        CHECK_EQUAL(expected, other);                            \
    }                                                            \
                                                                 \
    TEST_FIXTURE(Fixture, NAME ## _netural_random) {             \
        TEST_DETAILS();                                          \
        VECT neut, other = rand_vec(), expected;                 \
        for (KEY i=0; i<BASIS::dimension; ++i)                   \
            expected[i] OP other[i];                             \
                                                                 \
        neut OP other;                                           \
        CHECK_EQUAL(expected, neut);                             \
    }

TEST_INPLACE_BINARY_NEUTRAL(test_inplace_addition, +=)

TEST_INPLACE_BINARY_NEUTRAL(test_inplace_subtraction, -=)

#undef TEST_INPLACE_BINARY_NEUTRAL

#define TEST_INPLACE_BINARY_OPERATOR(NAME, OP)                   \
    TEST_FIXTURE(Fixture, NAME ## _random_random) {              \
        TEST_DETAILS();                                          \
        VECT lhs = rand_vec(), rhs = rand_vec(), expected(lhs);  \
        for (KEY i=0; i < BASIS::dimension; ++i)                 \
            expected[i] OP rhs[i];                               \
                                                                 \
        lhs OP rhs;                                              \
        CHECK_EQUAL(expected, lhs);                              \
    }

TEST_INPLACE_BINARY_OPERATOR(test_inplace_addition, +=)

TEST_INPLACE_BINARY_OPERATOR(test_inplace_subtraction, -=)


#undef TEST_INPLACE_BINARY_OPERATOR

TEST_FIXTURE(Fixture, test_inplace_scalar_mul_zero) {
    TEST_DETAILS();

    VECT vect(rand_vec()), expected;

    vect *= S(0);

    CHECK_EQUAL(expected, vect);
}

TEST_FIXTURE(Fixture, test_inplace_scalar_mul_1_ident) {
    TEST_DETAILS();
    VECT vect(rand_vec()), expected(vect);

    vect *= S(1);

    CHECK_EQUAL(expected, vect);
}

TEST_FIXTURE(Fixture, test_inplace_scalar_mul_mone) {
    TEST_DETAILS();
    VECT vect(rand_vec()), expected(-vect);
    vect *= S(-1);

    CHECK_EQUAL(expected, vect);
}


TEST_FIXTURE(Fixture, test_inplace_scalar_mul_random) {
    TEST_DETAILS();

    S sca (rand_scalar());
    VECT vect(rand_vec()), expected;

    for (KEY i=0; i < BASIS::dimension; ++i) {
        expected[i] = vect[i] * sca;
    }

    vect *= sca;
    CHECK_EQUAL(expected, vect);
}


TEST_FIXTURE(Fixture, test_inplace_rat_div_one_ident) {
    TEST_DETAILS();
    VECT vect(rand_vec()), expected(vect);

    vect /= Q(1);
    CHECK_EQUAL(expected, vect);
}

TEST_FIXTURE(Fixture, test_inplace_rat_div_mone) {
    TEST_DETAILS();
    VECT vect(rand_vec()), expected(-vect);

    vect /= Q(-1);
    CHECK_EQUAL(expected, vect);
}


TEST_FIXTURE(Fixture, test_inplace_rat_div_random) {
    TEST_DETAILS();

    Q rat (rand_scalar(S(1), S(5)));
    VECT vect(rand_vec()), expected;

    for (KEY i=0; i < BASIS::dimension; ++i) {
        expected[i] = vect[i] / rat;
    }

    vect /= rat;
    CHECK_EQUAL(expected, vect);
}

TEST_FIXTURE(Fixture, test_scalar_mul_zero) {
    TEST_DETAILS();

    VECT vect(rand_vec()), expected;

    CHECK_EQUAL(expected, vect * S(0));
}

TEST_FIXTURE(Fixture, test_scalar_mul_1_ident) {
    TEST_DETAILS();
    VECT vect(rand_vec()), expected(vect);

    CHECK_EQUAL(expected, vect * S(1));
}

TEST_FIXTURE(Fixture, test_scalar_mul_mone) {
    TEST_DETAILS();
    VECT vect(rand_vec()), expected(-vect);

    CHECK_EQUAL(expected, vect * S(-1));
}


TEST_FIXTURE(Fixture, test_scalar_mul_random) {
    TEST_DETAILS();

    S sca (rand_scalar());
    VECT vect(rand_vec()), expected;

    for (KEY i=0; i < BASIS::dimension; ++i) {
        expected[i] = vect[i] * sca;
    }

    CHECK_EQUAL(expected, vect * sca);
}


TEST_FIXTURE(Fixture, test_rat_div_one_ident) {
    TEST_DETAILS();
    VECT vect(rand_vec()), expected(vect);


    CHECK_EQUAL(expected, vect / Q(1));
}

TEST_FIXTURE(Fixture, test_rat_div_mone) {
    TEST_DETAILS();
    VECT vect(rand_vec()), expected(-vect);

    CHECK_EQUAL(expected, vect/ Q(-1));
}


TEST_FIXTURE(Fixture, test_rat_div_random) {
    TEST_DETAILS();

    Q rat (rand_scalar(S(1), S(5)));
    VECT vect(rand_vec()), expected;

    for (KEY i=0; i < BASIS::dimension; ++i) {
        expected[i] = vect[i] / rat;
    }

    CHECK_EQUAL(expected, vect / rat);
}

#define TEST_BINARY_OPERATOR_FUNC(NAME, OP, FUNC)                \
    TEST_FIXTURE(Fixture, NAME ## _random_random) {              \
        TEST_DETAILS();                                          \
        VECT lhs = rand_vec(), rhs = rand_vec(), expected;       \
        for (KEY i=0; i < BASIS::dimension; ++i)                 \
            expected[i] = FUNC (lhs[i], rhs[i]);                 \
                                                                 \
        CHECK_EQUAL(expected, lhs OP rhs);                       \
    }

TEST_BINARY_OPERATOR_FUNC(test_coordmin, &, std::min)

TEST_BINARY_OPERATOR_FUNC(test_coordmax, |, std::max)

#undef TEST_BINARY_OPERATOR_FUNC


#define TEST_INPLACE_BINARY_OPERATOR_FUNC(NAME, OP, FUNC)        \
    TEST_FIXTURE(Fixture, NAME ## _random_neutral) {             \
        TEST_DETAILS();                                          \
        VECT lhs(rand_vec()), rhs, expected;                     \
        for (KEY i=0; i<BASIS::dimension; ++i) {                 \
            std::pair<const KEY, S> p(i, FUNC(lhs[i], S(0)));    \
            expected.insert(p);                                  \
        }                                                        \
                                                                 \
        lhs OP rhs;                                              \
        CHECK_EQUAL(expected, lhs);                              \
    }                                                            \
                                                                 \
    TEST_FIXTURE(Fixture, NAME ## _neutral_random) {             \
        TEST_DETAILS();                                          \
        VECT lhs, rhs(rand_vec()), expected;                     \
        for (KEY i=0; i<BASIS::dimension; ++i) {                 \
            std::pair<const KEY, S> p(i, FUNC(S(0), rhs[i]));    \
            expected.insert(p);                                  \
        }                                                        \
                                                                 \
        lhs OP rhs;                                              \
        CHECK_EQUAL(expected, lhs);                              \
    }                                                            \
                                                                 \
    TEST_FIXTURE(Fixture, NAME ## _random_random) {              \
        TEST_DETAILS();                                          \
        VECT lhs = rand_vec(), rhs = rand_vec(), expected;       \
        for (KEY i=0; i < BASIS::dimension; ++i)  {              \
            std::pair<const KEY, S> p(i, FUNC(lhs[i], rhs[i]));  \
            expected.insert(p);                                  \
        }                                                        \
                                                                 \
        lhs OP rhs;                                              \
        CHECK_EQUAL(expected, lhs);                              \
    }

TEST_INPLACE_BINARY_OPERATOR_FUNC(test_inplace_coordmin, &=, std::min)

TEST_INPLACE_BINARY_OPERATOR_FUNC(test_inplace_coordmax, |=, std::max)

#undef TEST_INPLACE_BINARY_OPERATOR_FUNC


#define TEST_FUSED_OP(NAME, OP1, OP2, STYPE) \
    TEST_FIXTURE(Fixture, test_ ## NAME ## _vec) { \
        TEST_DETAILS();                             \
        VECT lhs = rand_vec(), rhs = rand_vec(), expected(lhs); \
        STYPE scal = (STYPE) rand_scalar();    \
        for (KEY i=0; i<BASIS::dimension; ++i)     \
            expected[i] OP1 (rhs[i] OP2 scal);     \
                                             \
        lhs. NAME(rhs, scal);                   \
        CHECK_EQUAL(expected, lhs);                         \
    }                                        \
                                             \
    TEST_FIXTURE(Fixture, test_ ## NAME ## _key) { \
        TEST_DETAILS();                      \
        VECT lhs = rand_vec(), expected(lhs);\
        KEY rhs = rand_key();                \
        STYPE scal = (STYPE) rand_scalar();    \
        expected[rhs] OP1 (S(1) OP2 scal);   \
                                             \
        lhs. NAME(rhs, scal);                \
        CHECK_EQUAL(expected, lhs);          \
    }

TEST_FUSED_OP(add_scal_prod, +=, *, S)

TEST_FUSED_OP(sub_scal_prod, -=, *, S)

TEST_FUSED_OP(add_scal_div, +=, /, Q)

TEST_FUSED_OP(sub_scal_div, -=, /, Q)

#undef TEST_FUSED_OP

TEST_FIXTURE (Fixture, test_unary_minus_neutral) {
    TEST_DETAILS();
    VECT neut, expected;

    CHECK_EQUAL(expected, -neut);
}

TEST_FIXTURE (Fixture, test_unary_minus_random) {
    TEST_DETAILS();
    VECT vect = rand_vec(), expected;
    for (KEY i = 0; i < BASIS::dimension; ++i)
        expected[i] = -vect[i];

    CHECK_EQUAL(expected, -vect);
}

TEST_FIXTURE(Fixture, test_unary_minus_plus_original_gives_netural) {
    TEST_DETAILS();
    VECT lhs, rhs(rand_vec()), expected;

    lhs = -rhs;

    CHECK_EQUAL(expected, lhs + rhs);
}

TEST_FIXTURE(Fixture, test_random_original_plus_unary_minus_gives_netural) {
    TEST_DETAILS();
    VECT lhs(rand_vec()), rhs(-lhs), expected;

    CHECK_EQUAL(expected, lhs + rhs);
}

#endif  // ifdef _VECTOR_TYPE
