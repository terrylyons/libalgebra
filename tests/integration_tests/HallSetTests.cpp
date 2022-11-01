// HallSetTests.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>

// the libalgebra framework
#include "alg_framework.h"
#include <libalgebra/libalgebra.h>

// Clang complains that operator<< is not declared for some reason,
// so forward declare here.
using alg::LET;
std::ostream& operator<<(std::ostream& os, std::pair<LET, LET> k);

// the unit test framework
#include "../common/time_and_details.h"
#include <UnitTest++.h>

#include "../common/compat.h"

std::ostream& operator<<(std::ostream& os, std::pair<LET, LET> k)
{
    return os << "(" + std::to_string(k.first) << ", " << std::to_string(k.second) << ")";
}

// validates the hall set and lie multiplication over it
SUITE(hallset)
{
    // DEPTH, ALPHABET SIZE, SCALAR TYPE
    typedef alg_framework<7, 5, Rational> SETUP75;// takes 143.6  seconds
    typedef alg_framework<6, 5, Rational> SETUP65;// takes   7.76 seconds
    typedef alg_framework<5, 7, Rational> SETUP57;// takes  10.89 seconds

    TEST_FIXTURE(SETUP65, hall_set_definition)
    {
        TEST_DETAILS();
        // the basis type
        typedef typename LIE::BASIS BASIS;
        // the basis instance
        const BASIS& basis = LIE::basis;
        const BASIS& hall_set = LIE::basis;
        // the free lie algebra is a vector space built over a hall basis
        // the hall basis is built from a collection of binary trees called a hall set
        // an instance of this basis type is described via an totally ordered sequence of KEYS of type LET starting at 1
        // typedef typename BASIS::KEY KEY;
        // where each key is connected to exactly two earlier keys, the key's parents
        // typedef typename BASIS::PARENT PARENTS;
        // each hall basis instance contains one (such tree defining) hall set
        using parents_t = typename BASIS::PARENT;

        // check the axioms of a hall set https://www.encyclopediaofmath.org/index.php/Hall_set

        // contains the god element
        CHECK_EQUAL(parents_t(0, 0), hall_set[0]);

        // contains A
        for (LET i = 1; i <= ALPHABET_SIZE; ++i)
            CHECK_EQUAL(parents_t(0, i), hall_set[i]);

        // order property
        for (typename std::vector<parents_t>::const_iterator it(hall_set.parents_begin());
             it != hall_set.parents_end(); ++it)
            if (*it != parents_t(0, 0))
                CHECK(it->first < it->second);

        // husband younger than father in-law property
        for (typename std::vector<parents_t>::const_iterator it(hall_set.parents_begin());
             it != hall_set.parents_end(); ++it)
            CHECK(it->first >= hall_set[it->second].first);

        // sanity tests

        // non-leaf strictly younger property (sanity test - follows from the above)
        for (LET i = 1 + ALPHABET_SIZE; i < hall_set.size(); ++i)
            CHECK((hall_set[i].second) < i);

        // test degree additivity and monotonicity
        // initial values
        for (LET i = 0; i < 1; ++i)
            CHECK_EQUAL(0U, basis.degree(i));
        for (LET i = 1; i < ALPHABET_SIZE + 1; ++i)
            CHECK_EQUAL(1U, basis.degree(i));
        // additivity
        for (LET i = 1; i < hall_set.size(); ++i)
            CHECK_EQUAL(basis.degree(hall_set[i].first) + basis.degree(hall_set[i].second), basis.degree(i));
        // monotonicity
        for (LET i = 1; i < hall_set.size(); ++i)
            CHECK(basis.degree(i - 1) <= basis.degree(i));

        // test grown_up enough
        CHECK(basis.degree(hall_set.size() - 1) >= DEPTH);

        // Check integrity of lie product and the reverse map on hall elements

        {
            typename std::vector<parents_t>::const_iterator begin = hall_set.parents_begin();
            typename std::map<parents_t, LET>::const_iterator it;

            // every hall set element is in the reverse map
            for (LET i = 1; i < hall_set.size(); ++i) {
                parents_t parents = *(begin + i);
                it = hall_set.find(parents);
                CHECK(it != hall_set.reverse_map_end());
                if (basis.degree(i) <= DEPTH) {
                    CHECK(LIE(i) == LIE(parents.first) * LIE(parents.second));
                }
                else {
                    CHECK(LIE() == LIE(parents.first) * LIE(parents.second));
                }
            }
        }

        // COMPUTATIONALLY HEAVY:
        // exhaustively confirm recursion for prod always terminates

        if (hall_set.size() * hall_set.size() < (1ULL << 31))
            for (LET i = 1; i < hall_set.size(); ++i)
                for (LET j = i + 1; j < hall_set.size(); ++j) {
                    // construct from LET type with constant one
                    LIE k = LIE(i) * LIE(j);// NON empty homogeneous

                    if (basis.degree(i) + basis.degree(j) > DEPTH)
                        CHECK(LIE() == k);
                    else {
                        CHECK(k.size() > 0);
                        for (typename LIE::iterator p(k.begin()); p != k.end(); ++p)
                            CHECK(basis.degree(iter::key<LIE>(p))
                                  == basis.degree(i) + basis.degree(j));
                    }

                    // check that items not in the hall basis are not in reversemap
                    if (i < hall_set[j].first) {
                        parents_t parents = parents_t(i, j);
                        typename std::map<parents_t, LET>::const_iterator it;
                        it = hall_set.find(parents);
                        CHECK(it == hall_set.reverse_map_end());
                    }
                }
    }
}
