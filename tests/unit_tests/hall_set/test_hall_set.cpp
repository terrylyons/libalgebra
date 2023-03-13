
//
// Created by sam on 07/10/2021.
//

#include <UnitTest++.h>
#include <libalgebra/libalgebra.h>
#include <utility>

struct hall_set_test_fixture {
    using hall_set_t = alg::hall_set<5>;

    using size_type     = typename hall_set_t::size_type;
    using degree_type   = typename hall_set_t::degree_type;
    using letter_type   = typename hall_set_t::letter_type;
    using key_type      = typename hall_set_t::key_type;
    using parent_type   = typename hall_set_t::parent_type;

    hall_set_t hall_set;

    hall_set_test_fixture() : hall_set{0}
    {}

    void grow_up(degree_type deg)
    {
        hall_set_t tmp(deg);
    }

    size_type start_of_degree(degree_type deg)
    {
        return hall_set_t::start_of_degree(deg);
    }

};



SUITE (hall_set_properties) {

    TEST_FIXTURE(hall_set_test_fixture, test_size_absolute) {
        // Out of order execution might cause this Hall set to already be larger than
        // degree 1, so we can only test that it is larger than each size.
        CHECK(5 <= hall_set.size());

        grow_up(2);
        CHECK(5 + 10 <= hall_set.size());

        grow_up(3);
        CHECK(5 + 10 + 40 <= hall_set.size());

        grow_up(4);
        CHECK(5 + 10 + 40 + 150 <= hall_set.size());
    }


    TEST_FIXTURE(hall_set_test_fixture, test_letters_size) {
        CHECK_EQUAL(5, hall_set.letters().size());
    }

    TEST_FIXTURE(hall_set_test_fixture, test_hall_set_degree_ranges) {
        grow_up(5);
        const auto& degree_ranges = hall_set.hall_set_degree_ranges();

        CHECK_EQUAL(size_type(0), degree_ranges[0].first);
        CHECK_EQUAL(size_type(1), degree_ranges[0].second);
        CHECK_EQUAL(size_type(1), degree_ranges[1].first);
        CHECK_EQUAL(size_type(6), degree_ranges[1].second);
        CHECK_EQUAL(size_type(6), degree_ranges[2].first);
        CHECK_EQUAL(size_type(16), degree_ranges[2].second);
        CHECK_EQUAL(size_type(16), degree_ranges[3].first);
        CHECK_EQUAL(size_type(56), degree_ranges[3].second);
        CHECK_EQUAL(size_type(56), degree_ranges[4].first);
        CHECK_EQUAL(size_type(206), degree_ranges[4].second);
        CHECK_EQUAL(size_type(206), degree_ranges[5].first);
        CHECK_EQUAL(size_type(830), degree_ranges[5].second);
    }

    TEST_FIXTURE(hall_set_test_fixture, test_letter_key_association) {
        for (const auto& let : hall_set.letters()) {
            CHECK(hall_set.letter(let));
            CHECK_EQUAL(let, hall_set.getletter(hall_set.keyofletter(let)));
        }
    }

    TEST_FIXTURE(hall_set_test_fixture, test_letter_lparent_rparent) {
        for (const auto& let: hall_set.letters()) {
            CHECK_EQUAL(letter_type(0), hall_set.lparent(let));
            CHECK_EQUAL(let, hall_set.rparent(let));
        }
    }

    TEST_FIXTURE(hall_set_test_fixture, test_lrparent_vs_lookup) {
        for (key_type k=1; k<hall_set.size(); ++k) {
            auto parents = hall_set[k];
            CHECK_EQUAL(parents.first, hall_set.lparent(k));
            CHECK_EQUAL(parents.second, hall_set.rparent(k));
        }
    }

    TEST_FIXTURE(hall_set_test_fixture, test_non_letter_letter_fail) {
        for (key_type k=6; k<hall_set.size(); ++k) {
            CHECK(!hall_set.letter(k));
        }
    }

    TEST_FIXTURE(hall_set_test_fixture, test_key_to_index) {
        for (key_type k=1; k<hall_set.size(); ++k) {
            CHECK_EQUAL(k-1, hall_set.key_to_index(k));
        }
    }

    TEST_FIXTURE(hall_set_test_fixture, test_key_to_index_roundtrip) {
        for (key_type k=1; k<hall_set.size(); ++k) {
            CHECK_EQUAL(k, hall_set.index_to_key(hall_set.key_to_index(k)));
        }
    }

    TEST_FIXTURE(hall_set_test_fixture, test_access_roundtrip) {
        for (key_type k=1; k<hall_set.size(); ++k) {
            CHECK_EQUAL(k, hall_set[hall_set[k]]);
        }
    }

    TEST_FIXTURE(hall_set_test_fixture, test_access_non_hall_member_fails) {
        parent_type p(2, 1); // not a member of the Hall set
        CHECK_THROW(hall_set[p], std::invalid_argument);
    }

    TEST_FIXTURE(hall_set_test_fixture, test_key_to_string_letters) {
        CHECK_EQUAL("1", hall_set.key2string(1));
        CHECK_EQUAL("2", hall_set.key2string(2));
        CHECK_EQUAL("3", hall_set.key2string(3));
        CHECK_EQUAL("4", hall_set.key2string(4));
        CHECK_EQUAL("5", hall_set.key2string(5));
    }

    TEST_FIXTURE(hall_set_test_fixture, test_key_to_string_deg2) {
        grow_up(2);
        CHECK_EQUAL("[1,2]", hall_set.key2string(6));
    }

    TEST_FIXTURE(hall_set_test_fixture, test_key_to_string_deg3) {
        grow_up(3);
        CHECK_EQUAL("[1,[1,2]]", hall_set.key2string(16));
    }


    TEST_FIXTURE(hall_set_test_fixture, test_begin_end_exhaustive) {
        grow_up(3); // make sure we have some elements

        key_type k1 = hall_set.begin();
        for (key_type k=1; k <= hall_set.size(); ++k) {
            CHECK_EQUAL(k, k1);
            k1 = hall_set.next_key(k1);
        }
        CHECK_EQUAL(key_type(0), k1);
    }



    TEST_FIXTURE(hall_set_test_fixture, test_hall_set_order)
    {
        grow_up(9);
        auto deg_1_range = hall_set.hall_set_degree_ranges()[1];
        for (key_type k=deg_1_range.first; k<deg_1_range.second; ++k) {
            // Degree 1 elements are letters.
            CHECK(hall_set.letter(k));
        }

        auto deg_2_range = hall_set.hall_set_degree_ranges()[2];
        for (key_type k=deg_2_range.first; k<deg_2_range.second; ++k) {
            // Product ab in H if and only if a and b are letters and a < b
            auto parents = hall_set[k];
            CHECK(hall_set.letter(parents.first));
            CHECK(hall_set.letter(parents.second));
            CHECK(parents.first < parents.second);
        }

        for (degree_type deg=3; deg<=hall_set.current_degree(); ++deg) {
            auto range = hall_set.hall_set_degree_ranges()[deg];
            for (key_type k=range.first; k<range.second; ++k) {
                auto parents = hall_set[k];
                auto rparents = hall_set[parents.second];
                auto a = parents.first;
                auto b = rparents.first;
                auto c = rparents.second;
                auto bc = parents.second;

                // Condition from Bourbaki is that
                // k in H if and only if k = a(bc) with bc in H and b <= a < bc and b < c
                CHECK(b <= a);
                CHECK(a < bc);
                CHECK(b < c);
            }
        }
    }


    TEST_FIXTURE(hall_set_test_fixture, test_hall_set_all_products_degree_2) {
        grow_up(2);

        const auto& letters = hall_set.letters();
        for (auto x : letters) {
            for (auto y : letters) {
                parent_type parents(x, y);
                if (x < y) {
                    auto xy = hall_set[parents];
                    auto reversed = hall_set[xy];
                    CHECK_EQUAL(x, reversed.first);
                    CHECK_EQUAL(y, reversed.second);
                } else {
                    CHECK_THROW(hall_set[parents], std::invalid_argument);
                }
            }
        }
    }

    TEST_FIXTURE(hall_set_test_fixture, test_hall_set_all_products_degree_3) {
        grow_up(3);

        const auto& letters = hall_set.letters();

        for (auto a : letters) {
            for (auto b : letters) {
                if (b <= a) {
                    for (auto c : letters) {
                        if (b < c) {
                            auto bc = hall_set[parent_type(b, c)];
                            if (a < bc) {
                                auto a_bc = hall_set[parent_type(a, bc)];
                                auto parents = hall_set[a_bc];
                                CHECK_EQUAL(a, parents.first);
                                CHECK_EQUAL(bc, parents.second);
                            } else {
                                CHECK_THROW(hall_set[parent_type(a, bc)], std::invalid_argument);
                            }
                        }
                    }
                }
            }
        }




    }




}
