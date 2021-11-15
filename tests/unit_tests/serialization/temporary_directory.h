//
// Created by sam on 26/10/2021.
//

#ifndef LIBALGEBRA_TEMPORARY_DIRECTORY_H
#define LIBALGEBRA_TEMPORARY_DIRECTORY_H

#include <boost/filesystem.hpp>

class temporary_directory : public boost::filesystem::path
{
    using path = boost::filesystem::path;

    temporary_directory();

public:
    ~temporary_directory();

    static temporary_directory& get_directory();
};

#endif//LIBALGEBRA_TEMPORARY_DIRECTORY_H
