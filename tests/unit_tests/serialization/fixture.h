//
// Created by sam on 26/10/2021.
//

#ifndef LIBALGEBRA_FIXTURE_H
#define LIBALGEBRA_FIXTURE_H

#include <fstream>
#include <iostream>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/filesystem.hpp>

struct fixture_base {

    using path = boost::filesystem::path;

    path archive;

    fixture_base();

    std::ofstream open_write() const;
    std::ifstream open_read() const;

    template<typename T>
    void write(const T& arg) const
    {
        auto file = open_write();
        boost::archive::text_oarchive oa{file};
        oa << arg;
    }

    template<typename T>
    T read() const
    {
        auto file = open_read();
        boost::archive::text_iarchive ia{file};
        T arg;
        ia >> arg;
        return arg;
    }
};

#endif//LIBALGEBRA_FIXTURE_H
