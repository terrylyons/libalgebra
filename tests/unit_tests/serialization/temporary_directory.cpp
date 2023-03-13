//
// Created by sam on 26/10/2021.
//

#include "temporary_directory.h"

temporary_directory::temporary_directory() : path(boost::filesystem::temp_directory_path() / boost::filesystem::unique_path())
{
    if (!boost::filesystem::create_directories(*this)) {
        throw std::runtime_error("Unable to create temporary directory");
    }
}

temporary_directory::~temporary_directory()
{
    boost::filesystem::remove_all(*this);
}

temporary_directory& temporary_directory::get_directory()
{
    static temporary_directory dir;
    return dir;
}
