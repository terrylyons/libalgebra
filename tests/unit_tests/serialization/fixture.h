//
// Created by sam on 26/10/2021.
//

#ifndef LIBALGEBRA_FIXTURE_H
#define LIBALGEBRA_FIXTURE_H

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/filesystem.hpp>

struct fixture_base {
    using path = boost::filesystem::path;

    path tmpdir;

    fixture_base() : tmpdir(boost::filesystem::temp_directory_path() / boost::filesystem::unique_path())
    {
        if (!boost::filesystem::create_directories(tmpdir)) {
            throw std::runtime_error("Unable to create temporary directory");
        }
    }

    ~fixture_base()
    {
        boost::filesystem::remove_all(tmpdir);
    }
};

#endif//LIBALGEBRA_FIXTURE_H
