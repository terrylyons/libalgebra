//
// Created by sam on 23/02/2021.
//

#if defined(_VECTOR_TYPE)

/*
 * Since our simple integer basis does not have degree, we are going to have
 * to use the tensor basis for our triangular binary transform methods.
 * We need to give the simple integer basis a multiplication operator in the
 * form of the following. We'll use tensor multiplication (concatenation) for
 * the other product.
 */


struct key_transform
{
    template <typename Vector, typename Key, typename Scalar>
    void operator()(Vector& result, const Key& k1, const Scalar& s1, const Key& k2, const Scalar& s2)
    {
        // diagonal multiply
        if (k1 == k2) {
            result[k1] += s1 * s2;
        }
    }
};

struct index_transform
{
    template <typename Scalar>
    void operator()(Scalar* result_ptr, const Scalar* lhs_ptr, const Scalar* rhs_ptr,
            const size_t lhs_target, const size_t rhs_target)
    {
        // Only if degree is the same
        // Note, for tensor basis the size of degree is unique for each n_letters/deg pair.
        if (lhs_target == rhs_target) {
            for (size_t i=0; i<lhs_target; ++i) {
                *(result_ptr++) += *(lhs_ptr++) * (*(rhs_ptr++));
            }
        }
    }
};



TEST_FIXTURE(Fixture, test_square_apply_binary_transform_key_only_empty_lhs) {
    TEST_DETAILS();

    VECT lhs, rhs(rand_vec()), expected, result;

    key_transform kt;
    lhs.square_buffered_apply_binary_transform(result, rhs, kt);

    CHECK_EQUAL(expected, result);
}

TEST_FIXTURE(Fixture, test_square_apply_binary_transform_key_index_empty_lhs) {
    TEST_DETAILS();

    VECT lhs, rhs(rand_vec()), expected, result;

    key_transform kt;
    index_transform it;
    lhs.square_buffered_apply_binary_transform(result, rhs, kt, it);

    CHECK_EQUAL(expected, result);
}

TEST_FIXTURE(Fixture, test_square_apply_binary_transform_key_only_empty_rhs) {
    TEST_DETAILS();

    VECT lhs(rand_vec()), rhs, expected, result;

    key_transform kt;
    lhs.square_buffered_apply_binary_transform(result, rhs, kt);

            CHECK_EQUAL(expected, result);
}

TEST_FIXTURE(Fixture, test_square_apply_binary_transform_key_index_empty_rhs) {
    TEST_DETAILS();

    VECT lhs(rand_vec()), rhs, expected, result;

    key_transform kt;
    index_transform it;
    lhs.square_buffered_apply_binary_transform(result, rhs, kt, it);

    CHECK_EQUAL(expected, result);
}

TEST_FIXTURE(Fixture, test_square_apply_binary_transform_key_only_random) {
    TEST_DETAILS();

    VECT lhs(rand_vec()), rhs(rand_vec()), expected, result;

    for (KEY i=0; i < BASIS::dimension; ++i) {
        expected[i] = lhs[i] * rhs[i];
    }

    key_transform kt;
    lhs.square_buffered_apply_binary_transform(result, rhs, kt);

    CHECK_EQUAL(expected, result);
}

TEST_FIXTURE(Fixture, test_square_apply_binary_transform_key_index_random) {
    TEST_DETAILS();

    VECT lhs(rand_vec()), rhs(rand_vec()), expected, result;

    for (KEY i=0; i < BASIS::dimension; ++i) {
        expected[i] = lhs[i] * rhs[i];
    }

    key_transform kt;
    index_transform it;
    lhs.square_buffered_apply_binary_transform(result, rhs, kt, it);

    CHECK_EQUAL(expected, result);
}

TEST_FIXTURE(Fixture, test_square_apply_binary_transform_key_only_unidim_different) {
    TEST_DETAILS();

    VECT lhs(KEY(0), rand_scalar()), rhs(KEY(1), rand_scalar()), result, expected;

    key_transform kt;
    lhs.square_buffered_apply_binary_transform(result, rhs, kt);

    CHECK_EQUAL(expected, result);
}

TEST_FIXTURE(Fixture, test_square_apply_binary_transform_key_only_unidim_same) {
    TEST_DETAILS();

    VECT lhs(KEY(1), rand_scalar()), rhs(KEY(1), rand_scalar()), result, expected;
    expected[KEY(1)] = lhs[KEY(1)] * rhs[KEY(1)];

    key_transform kt;
    lhs.square_buffered_apply_binary_transform(result, rhs, kt);

    CHECK_EQUAL(expected, result);
}


TEST_FIXTURE(Fixture, test_square_apply_binary_transform_key_index_unidim_same) {
    TEST_DETAILS();

    VECT lhs(KEY(1), rand_scalar()), rhs(KEY(1), rand_scalar()), result, expected;
    expected[KEY(1)] = lhs[KEY(1)] * rhs[KEY(1)];

    key_transform kt;
    index_transform it;
    lhs.square_buffered_apply_binary_transform(result, rhs, kt, it);

    CHECK_EQUAL(expected, result);
}


TEST_FIXTURE(Fixture, test_square_apply_binary_transform_key_only_existing_data) {
    TEST_DETAILS();

    VECT lhs(rand_vec()), rhs(rand_vec()), result(rand_vec()), expected;
    for (KEY i=0; i < BASIS::dimension; ++i) {
        expected[i] = result[i] + (lhs[i] * rhs[i]);
    }

    key_transform kt;
    lhs.square_buffered_apply_binary_transform(result, rhs, kt);

    CHECK_EQUAL(expected, result);
}


TEST_FIXTURE(Fixture, test_square_apply_binary_transform_key_index_existing_data) {
    TEST_DETAILS();

    VECT lhs(rand_vec()), rhs(rand_vec()), result(rand_vec()), expected;
    for (KEY i=0; i < BASIS::dimension; ++i) {
        expected[i] = result[i] + (lhs[i] * rhs[i]);
    }

    key_transform kt;
    index_transform it;
    lhs.square_buffered_apply_binary_transform(result, rhs, kt, it);

    CHECK_EQUAL(expected, result);
}


struct ident
{
    template <typename Scalar>
    Scalar operator()(Scalar v) { return v; }
};


TEST_FIXTURE(Fixture, test_triangular_apply_binary_transform_key_only_empty_lhs) {
    TEST_DETAILS();

    TVECT lhs, rhs(rand_tvec(2)), expected, result;

    ft_key_operator<ident> kt {ident()};
    lhs.triangular_buffered_apply_binary_transform(result, rhs, kt, TBASIS::MAX_DEGREE);

    CHECK_EQUAL(expected.size(), result.size());
    CHECK_EQUAL(expected, result);
}

TEST_FIXTURE(Fixture, test_triangular_apply_binary_transform_key_index_empty_lhs) {
    TEST_DETAILS();

    TVECT lhs, rhs(rand_tvec(2)), expected, result;

    ft_key_operator<ident> kt {ident()};
    ft_index_operator<ident> it{ident()};
    lhs.triangular_buffered_apply_binary_transform(result, rhs, kt, it, TBASIS::MAX_DEGREE);

    CHECK_EQUAL(expected.size(), result.size());
    CHECK_EQUAL(expected, result);
}

TEST_FIXTURE(Fixture, test_triangular_apply_binary_transform_key_only_empty_rhs) {
    TEST_DETAILS();

    TVECT lhs(rand_tvec(2)), rhs, expected, result;

    ft_key_operator<ident> kt{ident()};
    lhs.triangular_buffered_apply_binary_transform(result, rhs, kt, TBASIS::MAX_DEGREE);

    CHECK_EQUAL(expected.size(), result.size());
    CHECK_EQUAL(expected, result);
}

TEST_FIXTURE(Fixture, test_triangular_apply_binary_transform_key_index_empty_rhs) {
    TEST_DETAILS();

    TVECT lhs(rand_tvec(2)), rhs, expected, result;

    ft_key_operator<ident> kt{ident()};
    ft_index_operator<ident> it{ident()};
    lhs.triangular_buffered_apply_binary_transform(result, rhs, kt, it,  TBASIS::MAX_DEGREE);

    CHECK_EQUAL(expected.size(), result.size());
    CHECK_EQUAL(expected, result);
}

TEST_FIXTURE(Fixture, test_triangular_apply_binary_transform_key_only_random) {
    TEST_DETAILS();

    TVECT lhs(rand_tvec(2)), rhs(rand_tvec(2)), expected, result;
    TKEY l, k, empty;
    while (result.basis.degree(k) <= 2) {
        while (result.basis.degree(l) <= 2) {
            expected[k * l] += lhs[k] * rhs[l];
            l = result.basis.nextkey(l);
        }
        l = empty;
        k = result.basis.nextkey(k);
    }

    ft_key_operator<ident> kt{ident()};
    lhs.triangular_buffered_apply_binary_transform(result, rhs, kt, TBASIS::MAX_DEGREE);

    CHECK_EQUAL(expected.size(), result.size());
    CHECK_EQUAL(expected, result);
}

TEST_FIXTURE(Fixture, test_triangular_apply_binary_transform_key_index_random) {
    TEST_DETAILS();

    TVECT lhs(rand_tvec(2)), rhs(rand_tvec(2)), expected, result;

    TKEY l, k, empty;
    while (result.basis.degree(k) <= 2) {
        while (result.basis.degree(l) <= 2) {
            expected[k * l] += lhs[k] * rhs[l];
            l = result.basis.nextkey(l);
        }
        l = empty;
        k = result.basis.nextkey(k);
    }

    ft_key_operator<ident> kt{ident()};
    ft_index_operator<ident> it{ident()};
    lhs.triangular_buffered_apply_binary_transform(result, rhs, kt, it, TBASIS::MAX_DEGREE);

    CHECK_EQUAL(expected.size(), result.size());
    CHECK_EQUAL(expected, result);
}

TEST_FIXTURE(Fixture, test_triangular_apply_binary_transform_key_only_unidim_different) {
    TEST_DETAILS();

    TVECT lhs(TKEY(LET(1)), rand_scalar()), rhs(TKEY(LET(2)), rand_scalar()), result, expected;

    expected[TKEY(LET(1)) * TKEY(LET(2))] = lhs[TKEY(LET(1))] * rhs[TKEY(LET(2))];

    ft_key_operator<ident> kt{ident()};
    lhs.triangular_buffered_apply_binary_transform(result, rhs, kt,  TBASIS::MAX_DEGREE);

    CHECK_EQUAL(expected.size(), result.size());
    CHECK_EQUAL(expected, result);
}

TEST_FIXTURE(Fixture, test_triangular_apply_binary_transform_key_only_unidim_same) {
    TEST_DETAILS();

    TVECT lhs(TKEY(LET(1)), rand_scalar()), rhs(TKEY(LET(1)), rand_scalar()), result, expected;
    expected[TKEY(LET(1))*TKEY(LET(1))] = lhs[TKEY(LET(1))] * rhs[TKEY(LET(1))];

    ft_key_operator<ident> kt{ident()};
    lhs.triangular_buffered_apply_binary_transform(result, rhs, kt, TBASIS::MAX_DEGREE);

    CHECK_EQUAL(expected.size(), result.size());
    CHECK_EQUAL(expected, result);
}

TEST_FIXTURE(Fixture, test_triangular_apply_binary_transform_key_index_unidim_different) {
    TEST_DETAILS();

    TVECT lhs(TKEY(LET(1)), rand_scalar()), rhs(TKEY(LET(2)), rand_scalar()), result, expected;
    expected[TKEY(LET(1))*TKEY(LET(2))] = lhs[TKEY(LET(1))] * rhs[TKEY(LET(2))];

    ft_key_operator<ident> kt{ident()};
    ft_index_operator<ident> it{ident()};
    lhs.triangular_buffered_apply_binary_transform(result, rhs, kt, it, TBASIS::MAX_DEGREE);

    CHECK_EQUAL(expected.size(), result.size());
    CHECK_EQUAL(expected, result);
}

TEST_FIXTURE(Fixture, test_triangular_apply_binary_transform_key_index_unidim_same) {
    TEST_DETAILS();

    TVECT lhs(TKEY(LET(1)), rand_scalar()), rhs(TKEY(LET(1)), rand_scalar()), result, expected;
    expected[TKEY(LET(1))*TKEY(LET(1))] = lhs[TKEY(LET(1))] * rhs[TKEY(LET(1))];

    ft_key_operator<ident> kt{ident()};
    ft_index_operator<ident> it{ident()};
    lhs.triangular_buffered_apply_binary_transform(result, rhs, kt, it, TBASIS::MAX_DEGREE);

    CHECK_EQUAL(expected.size(), result.size());
    CHECK_EQUAL(expected, result);
}


TEST_FIXTURE(Fixture, test_triangular_apply_binary_transform_key_only_existing_data) {
    TEST_DETAILS();

    TVECT lhs(rand_tvec(2)), rhs(rand_tvec(2)), result(rand_tvec(2)), expected(result);
    TKEY l, k, empty;
    while (result.basis.degree(k) <= 2) {
        while (result.basis.degree(l) <= 2) {
            expected[k * l] += (lhs[k] * rhs[l]);
            l = result.basis.nextkey(l);
        }
        l = empty;
        k = result.basis.nextkey(k);
    }

    ft_key_operator<ident> kt{ident()};
    lhs.triangular_buffered_apply_binary_transform(result, rhs, kt, TBASIS::MAX_DEGREE);


    CHECK_EQUAL(expected.size(), result.size());
    CHECK_EQUAL(expected, result);
}


TEST_FIXTURE(Fixture, test_triangular_apply_binary_transform_key_index_existing_data) {
    TEST_DETAILS();

    TVECT lhs(rand_tvec(2)), rhs(rand_tvec(2)), result(rand_tvec(2)), expected(result);
    TKEY l, k, empty;
    while (result.basis.degree(k) <= 2) {
        while (result.basis.degree(l) <= 2) {
            expected[k * l] +=  (lhs[k] * rhs[l]);
            l = result.basis.nextkey(l);
        }
        l = empty;
        k = result.basis.nextkey(k);
    }

    ft_key_operator<ident> kt{ident()};
    ft_index_operator<ident> it{ident()};
    lhs.triangular_buffered_apply_binary_transform(result, rhs, kt, it, TBASIS::MAX_DEGREE);

    CHECK_EQUAL(expected.size(), result.size());
    CHECK_EQUAL(expected, result);
}


#endif