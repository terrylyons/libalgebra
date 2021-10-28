//
// Created by sam on 28/10/2021.
//

#include <UnitTest++/UnitTest++.h>
#include <libalgebra/hall_set.h>

#include <fstream>
#include <iostream>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/filesystem.hpp>

struct fixture {

    using hall_set_t = alg::hall_set<17>;

    alg::hall_set<17> hall_set;

    using size_type = typename hall_set_t::size_type;
    using degree_type = typename hall_set_t::degree_type;
    using letter_type = typename hall_set_t::letter_type;
    using key_type = typename hall_set_t::key_type;
    using parent_type = typename hall_set_t::parent_type;

    fixture() : hall_set(0)
    {
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
            oa << tmp;
        }
    }
};

SUITE(hall_set_serialization)
{

    TEST_FIXTURE(fixture, test_letters_size)
    {
        CHECK_EQUAL(size_type(17), hall_set.letters().size());
    }

    TEST_FIXTURE(fixture, test_current_degree)
    {
        CHECK_EQUAL(degree_type(2), hall_set.current_degree());
    }

    TEST_FIXTURE(fixture, test_size_of_data)
    {
        CHECK_EQUAL(size_type(0), hall_set.size());
    }


}