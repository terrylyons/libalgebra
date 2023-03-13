//
// Created by sam on 28/10/2021.
//

#include <UnitTest++.h>
#include <libalgebra/hall_set.h>

#include <fstream>
#include <iostream>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/filesystem.hpp>

SUITE(hall_set_serialization)
{

    struct fixture {

        using hall_set_t = alg::hall_set<17>;

        hall_set_t hall_set;
        static bool loaded;

        fixture() : hall_set(1)
        {
        }

        using size_type = typename hall_set_t::size_type;
        using degree_type = typename hall_set_t::degree_type;
        using letter_type = typename hall_set_t::letter_type;
        using key_type = typename hall_set_t::key_type;
        using parent_type = typename hall_set_t::parent_type;

        void load()
        {
            if (loaded) return;
            loaded = true;
            boost::filesystem::path path{"hall_set-17-2.txt"};
            if (exists(path)) {
                std::ifstream is(path.c_str());
                boost::archive::text_iarchive ia(is);
                ia >> hall_set;
            }
            else {
                alg::hall_set<17> tmp(2);
                std::ofstream os(path.c_str());
                boost::archive::text_oarchive oa(os);
                oa << hall_set;
            }
        }
    };

    bool fixture::loaded = false;

    TEST_FIXTURE(fixture, test_letters_size)
    {
        load();
        CHECK_EQUAL(size_type(17), hall_set.letters().size());
    }

    TEST_FIXTURE(fixture, test_current_degree)
    {
        load();
        CHECK_EQUAL(degree_type(2), hall_set.current_degree());
    }

    TEST_FIXTURE(fixture, test_size_of_data)
    {
        load();
        CHECK_EQUAL(size_type(153), hall_set.size());
    }
}
