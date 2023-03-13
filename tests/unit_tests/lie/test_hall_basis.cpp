//
// Created by sam on 02/02/2021.
//

#include <UnitTest++.h>

#include <libalgebra/libalgebra.h>
#include <libalgebra/alg_types.h>

#include "../../common/time_and_details.h"


using alg::DEG;
using alg::LET;



SUITE(hall_basis) {


    const unsigned test_alphabet_size = 5;
    const unsigned test_max_degree = 4;
    typedef double TestScalarType;
    typedef double TestRationalType;

    typedef alg::coefficients::coefficient_field<double> coeff_t;

    typedef alg::lie<coeff_t, test_alphabet_size, test_max_degree> LIE;
    typedef alg::lie_basis<test_alphabet_size, test_max_degree> LIEBASIS;

    class TestBasis : public alg::hall_basis<test_alphabet_size, test_max_degree>
    {
    public:
        typedef double RATIONAL;
        typedef double SCALAR;


        typedef typename alg::hall_basis<test_alphabet_size, test_max_degree>::KEY KEY;
        const std::vector<std::pair<size_t, size_t> >& get_hall_set_degree_ranges() const
        {
            return hall_set_degree_ranges();
        }

        const std::vector<LET>& get_letters() const
        {
            return letters();
        }

        const std::vector<KEY>& get_ltk() const
        {
            return l2k();
        }

    };

    typedef typename TestBasis::KEY KEY;


    TEST(test_basis_created_correct_num_letters) {
        TEST_DETAILS();
        TestBasis basis;

        CHECK_EQUAL(test_alphabet_size, basis.get_letters().size());
    }

    TEST(test_basis_created_ltk_size) {
        TEST_DETAILS();
        TestBasis basis;

        CHECK_EQUAL(test_alphabet_size, basis.get_ltk().size());
    }

/*
    TEST(test_growup_functionality) {
        TEST_DETAILS();
        TestBasis basis;

        CHECK_EQUAL(2, basis.get_hall_set_degree_ranges().size());

        basis.growup(2);

        CHECK_EQUAL(3, basis.get_hall_set_degree_ranges().size());
    }
    */

    TEST(test_degree) {
        TEST_DETAILS();
        TestBasis basis;


        CHECK_EQUAL(0, basis.degree(0));
        for (LET i=1; i<=test_alphabet_size; ++i) {
            CHECK_EQUAL(1, basis.degree(i));
        }
    }


    TEST(test_key_of_letter) {
        TEST_DETAILS();
        TestBasis basis;

        for (LET i=1; i<=test_alphabet_size; ++i) {
            CHECK_EQUAL(KEY(i), basis.keyofletter(i));
        }

    }

    TEST(test_lparent_of_key_base_letters) {
        TEST_DETAILS();
        TestBasis basis;


        for (LET i=1; i<=test_alphabet_size; ++i) {
            CHECK_EQUAL(0, basis.lparent(i));
        }

    }

    TEST(test_rparent_of_key_base_letters) {
        TEST_DETAILS();
        TestBasis basis;

        for (LET i=1; i<=test_alphabet_size; ++i) {
            CHECK_EQUAL(i, basis.rparent(i));
        }

    }

    TEST(test_lparent_of_key_higher_key) {
        TEST_DETAILS();
        TestBasis basis;

        KEY key = test_alphabet_size + 1;
        CHECK_EQUAL(1, basis.lparent(key));

    }

    TEST(test_rparent_of_key_higher_key) {
        TEST_DETAILS();
        TestBasis basis;

        KEY key = test_alphabet_size + 1;
        CHECK_EQUAL(2, basis.rparent(key));
    }

    TEST(test_is_letter_for_letters) {
        TEST_DETAILS();
        TestBasis basis;

        for (LET i=1; i<=test_alphabet_size; ++i) {
            CHECK(basis.letter(i));
        }
    }

    TEST(test_is_letter_non_letter) {
        TEST_DETAILS();
        TestBasis basis;

        KEY key = test_alphabet_size + 1;
        CHECK(!basis.letter(key));
    }

    TEST (test_begin_end_exhaustive)
    {
        TestBasis hall_set;
        KEY k1 = hall_set.begin();
        for (KEY k = 1; k<=hall_set.size(); ++k) {
            CHECK_EQUAL(k, k1);
            k1 = hall_set.nextkey(k1);
        }
        CHECK_EQUAL(KEY(0), k1);
    }

}



SUITE(hall_basis_size) {

#define TEST_HALL_SET_SIZE(WIDTH, DEPTH, EXPECTED)                                  \
    TEST(test_hall_set_size_ ## WIDTH ## _ ## DEPTH) {                              \
        TEST_DETAILS();                                                             \
        CHECK_EQUAL(EXPECTED, alg::dtl::hall_set_size(WIDTH, DEPTH));               \
    }

    TEST_HALL_SET_SIZE(2, 0, 0)
    TEST_HALL_SET_SIZE(2, 1, 2)
    TEST_HALL_SET_SIZE(2, 2, 3)
    TEST_HALL_SET_SIZE(2, 3, 5)
    TEST_HALL_SET_SIZE(2, 4, 8)
    TEST_HALL_SET_SIZE(2, 5, 14)
    TEST_HALL_SET_SIZE(2, 6, 23)
    TEST_HALL_SET_SIZE(2, 7, 41)
    TEST_HALL_SET_SIZE(2, 8, 71)
    TEST_HALL_SET_SIZE(2, 9, 127)


    TEST_HALL_SET_SIZE(3, 0, 0)
    TEST_HALL_SET_SIZE(3, 1, 3)
    TEST_HALL_SET_SIZE(3, 2, 6)
    TEST_HALL_SET_SIZE(3, 3, 14)
    TEST_HALL_SET_SIZE(3, 4, 32)
    TEST_HALL_SET_SIZE(3, 5, 80)
    TEST_HALL_SET_SIZE(3, 6, 196)
    TEST_HALL_SET_SIZE(3, 7, 508)
    TEST_HALL_SET_SIZE(3, 8, 1318)
    TEST_HALL_SET_SIZE(3, 9, 3502)


    TEST_HALL_SET_SIZE(4, 0, 0)
    TEST_HALL_SET_SIZE(4, 1, 4)
    TEST_HALL_SET_SIZE(4, 2, 10)
    TEST_HALL_SET_SIZE(4, 3, 30)
    TEST_HALL_SET_SIZE(4, 4, 90)
    TEST_HALL_SET_SIZE(4, 5, 294)
    TEST_HALL_SET_SIZE(4, 6, 964)
    TEST_HALL_SET_SIZE(4, 7, 3304)
    TEST_HALL_SET_SIZE(4, 8, 11464)
    TEST_HALL_SET_SIZE(4, 9, 40584)


    TEST_HALL_SET_SIZE(5, 0, 0)
    TEST_HALL_SET_SIZE(5, 1, 5)
    TEST_HALL_SET_SIZE(5, 2, 15)
    TEST_HALL_SET_SIZE(5, 3, 55)
    TEST_HALL_SET_SIZE(5, 4, 205)
    TEST_HALL_SET_SIZE(5, 5, 829)
    TEST_HALL_SET_SIZE(5, 6, 3409)
    TEST_HALL_SET_SIZE(5, 7, 14569)
    TEST_HALL_SET_SIZE(5, 8, 63319)
    TEST_HALL_SET_SIZE(5, 9, 280319)

#undef TEST_HALL_SET_SIZE

}
