//
// Created by sam on 26/10/2021.
//

#include "fixture.h"
#include "temporary_directory.h"

fixture_base::fixture_base() : archive(temporary_directory::get_directory() / boost::filesystem::unique_path())
{}

std::ifstream fixture_base::open_read() const
{
    return std::ifstream(archive.c_str());
}

std::ofstream fixture_base::open_write() const
{
    return std::ofstream(archive.c_str());
}