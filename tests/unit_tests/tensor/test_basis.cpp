//
// Created by sam on 02/02/2021.
//

#include <UnitTest++.h>

#include <libalgebra/coefficients.h>
#include <libalgebra/libalgebra.h>

#include "../../common/time_and_details.h"

using alg::DEG;
using alg::LET;

template<typename S, typename R, DEG W, DEG D>
struct BasisTool {
    typedef alg::coefficients::coefficient_ring<S, R> Coefficients;
    typedef alg::tensor_basis<W, D> TBASIS;
    typedef alg::free_tensor_basis<W, D> FTBASIS;
    typedef alg::shuffle_tensor_basis<W, D> STBASIS;
    typedef typename TBASIS::KEY KEY;

    static const DEG n_letters = W;
    static const DEG max_depth = D;

    KEY make_key(const LET* letters, std::size_t N)
    {
        KEY k;
        for (size_t i = 0; i < N; ++i)
            k.push_back(letters[i]);
        return k;
    }
};

struct Framework55Double : public BasisTool<double, double, 5, 5> {
    typedef BasisTool<double, double, 5, 5> TOOL;
    typedef typename TOOL::TBASIS TBASIS;
    typedef typename TOOL::FTBASIS FTBASIS;
    typedef typename TOOL::STBASIS STBASIS;
    typedef typename TOOL::KEY KEY;
};

SUITE(tensor_basis)
{

    TEST_FIXTURE(Framework55Double, test_degree)
    {
        TEST_DETAILS();

        KEY k;
        TBASIS basis;

        CHECK_EQUAL(0, basis.degree(k));

        for (DEG d = 1; d <= 5; ++d) {
            k.push_back(LET(1));
            CHECK_EQUAL(d, basis.degree(k));
        }
    }

    TEST_FIXTURE(Framework55Double, test_begin_key_empty)
    {
        TEST_DETAILS();
        TBASIS basis;

        CHECK_EQUAL(KEY(), basis.begin());
    }

    /*
    TEST_FIXTURE(Framework55Double, test_last_key) {
        TEST_DETAILS()

        TBASIS basis;

        std::array<LET, 5> letters  {4, 4, 4, 4, 4};
        CHECK_EQUAL(make_key(letters), basis.end());
    }
    */

    TEST_FIXTURE(Framework55Double, test_next_key)
    {
        TEST_DETAILS()

        TBASIS basis;

        LET in[] = {1, 2, 3};
        KEY k = make_key(in, 3);

        LET expected[] = {1, 2, 4};
        CHECK_EQUAL(make_key(expected, 3), basis.nextkey(k));
    }

    TEST_FIXTURE(Framework55Double, test_key_to_index_empty)
    {
        KEY k;

        FTBASIS basis;
        CHECK_EQUAL(0, basis.key_to_index(k));
    }

    TEST_FIXTURE(Framework55Double, test_key_to_index_letters)
    {
        KEY k;

        FTBASIS basis;

        for (LET i = 1; i <= n_letters; ++i) {
            CHECK_EQUAL(i, basis.key_to_index(KEY(i)));
        }
    }

    TEST_FIXTURE(Framework55Double, test_key_to_index_deg_2)
    {
        KEY k(KEY(LET(1)) * KEY(LET(1)));

        FTBASIS basis;

        size_t idx = n_letters;
        while (basis.degree(k) == 2) {
            CHECK_EQUAL(++idx, basis.key_to_index(k));
            k = basis.nextkey(k);
        }
    }
}

SUITE(free_tensor_basis)
{

    TEST_FIXTURE(Framework55Double, test_prod_empty_keys)
    {
        TEST_DETAILS()

        FTBASIS basis;
        KEY k1, k2, expected;

        CHECK_EQUAL(expected, basis.prod(k1, k2));
    }

    TEST_FIXTURE(Framework55Double, test_prod_empty_deg1)
    {
        TEST_DETAILS()

        FTBASIS basis;

        KEY k1;
        KEY expected = KEY(LET(1));

        CHECK_EQUAL(expected, basis.prod(KEY(LET(1)), k1));
        CHECK_EQUAL(expected, basis.prod(k1, KEY(LET(1))));
    }

    TEST_FIXTURE(Framework55Double, test_prod_deg1_deg1)
    {
        TEST_DETAILS()

        FTBASIS basis;

        LET expected_let[] = {1, 2};
        KEY expected = make_key(expected_let, 2);
        KEY result = basis.prod(KEY(LET(1)), KEY(LET(2)));

        CHECK_EQUAL((DEG)2, basis.degree(result));
        CHECK_EQUAL(expected, result);
    }

    TEST_FIXTURE(Framework55Double, test_prod_deg1_deg2)
    {
        TEST_DETAILS()

        FTBASIS basis;

        LET rhs_let[] = {2, 3};
        LET expected_let[] = {1, 2, 3};
        KEY expected = make_key(expected_let, 3);
        KEY result = basis.prod(KEY(LET(1)), make_key(rhs_let, 2));

        CHECK_EQUAL((DEG)3, basis.degree(result));
        CHECK_EQUAL(expected, result);
    }

    TEST_FIXTURE(Framework55Double, nextkey_TENSOR)
    {
        FTBASIS basis;

        //// enumerate the words manually into a vector and compare
        auto width = FTBASIS::NO_LETTERS;
        auto max_degree = FTBASIS::MAX_DEGREE;
        KEY k = basis.begin();
        std::vector<KEY> keys(1);

        CHECK_EQUAL(k, keys.back());
        k = basis.nextkey(k);

        size_t j = 0, j_max = keys.size();
        for (size_t depth = 1; depth <= max_degree; ++depth) {
            for (; j != j_max; ++j) {
                // join all letters i to each word j of one degree less
                for (LET i = 1; i <= width; ++i) {
                    // KEY(LET arg) requires 0 < arg <= NO_LETTERS
                    keys.push_back(keys[j] * KEY(i));
                    CHECK_EQUAL(k, keys.back());
                    k = basis.nextkey(k);
                    //// show words in order placed on vector
                }
            }
            j_max = keys.size();
        }

        CHECK_EQUAL(k, basis.end());
    }
}
