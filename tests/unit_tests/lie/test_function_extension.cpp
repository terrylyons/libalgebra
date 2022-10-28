//
// Created by sam on 25/10/2021.
//

#include <UnitTest++.h>
#include <libalgebra/libalgebra.h>

#include <libalgebra/alg_types.h>

struct Fixture : alg_types<5, 5, DPReal> {
    using base = alg_types<5, 5, DPReal>;
    using hall_basis_t = typename base::LIE::BASIS;
    using key_type = typename hall_basis_t::KEY;

    hall_basis_t hall_basis;

    static LIE multiindex_letter(LET l)
    {
        const auto& hb = LIE::basis;
        return LIE(hb.keyofletter(l));
    }

    static LIE multiindex_bino(const LIE& lhs, const LIE& rhs)
    {
        return lhs + rhs;
    }

    using multiindex_func_t = hall_basis_t::extended_function<
            decltype(&multiindex_letter),
            decltype(&multiindex_bino),
            alg::no_caching_tag>;

    Fixture() : hall_basis(), multiindex(hall_basis, multiindex_letter, multiindex_bino)
    {}

    multiindex_func_t multiindex;
};



SUITE(hall_extended_function) {

    TEST_FIXTURE(Fixture, test_extended_function_letters)
    {
        for (LET i=1; i<=5; ++i) {
            LIE expected(hall_basis.keyofletter(i));
            CHECK_EQUAL(expected, multiindex(key_type(i)));
        }
    }


    TEST_FIXTURE(Fixture, test_extended_function_deg_2)
    {
        key_type k = 6; // [1, 2]

        LIE expected;
        expected.add_scal_prod(key_type(1), 1.0);
        expected.add_scal_prod(key_type(2), 1.0);

        CHECK_EQUAL(expected, multiindex(k));
    }

    TEST_FIXTURE(Fixture, test_extended_function_deg_3)
    {
        auto k = static_cast<key_type>(hall_basis_t::start_of_degree(3) + 1); // [1, [1, 2]]
        LIE expected;
        expected.add_scal_prod(key_type(1), 2.0);
        expected.add_scal_prod(key_type(2), 1.0);

        CHECK_EQUAL(expected, multiindex(k));
    }



}
