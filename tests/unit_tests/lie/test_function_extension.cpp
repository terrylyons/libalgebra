//
// Created by sam on 25/10/2021.
//

#include <UnitTest++/UnitTest++.h>
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
