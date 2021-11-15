#include "memfile.h"

char* memfile::begin() const
{
	return b;
}

char* memfile::end() const
{
	return e;
}

const char* memfile::cbegin() const
{
	return b;
}

const char* memfile::cend() const
{
	return e;
}

size_t memfile::size() const
{
	return e - b;
}

memfile::memfile(boost::filesystem::path filename, size_t numberOfBytes /*= 0*/) : p(std::move(filename)), readonly(false)
{
	// populate file if it does not already exist
	//boost::filesystem::path p(filename.c_str());
	int count = 10;
	if (boost::filesystem::exists(p))
		read_only(true);
	else
	while (!boost::filesystem::exists(p) && (count-- > 0)) {
		// create file
		{
			// simply create the file
			if (0 == count) throw;// multiple attempts to create file failed
			boost::filesystem::ofstream f;
			f.open(p);
		}// file closed

		// only size the file on first use, then fixed
		boost::filesystem::resize_file(p, numberOfBytes);//62496
	}

	// open the file for read and write
    boost::iostreams::mapped_file_params params;
	params.path = p.string();
	params.flags = (!read_only()) ? boost::iostreams::mapped_file_sink::readwrite : boost::iostreams::mapped_file_sink::priv;
	file.open(params);
	if (!(file.is_open()))
		throw;
	// in readonly mode file is private - you can write to the memory but it will not be copied back to the file
	b = file.begin();
	e = file.end();
	//numberOfBytes = end() - begin();
}

memfile::~memfile()
{
	file.close();
}

bool memfile::read_only(bool val)
{
	return readonly = readonly || val;
}

bool memfile::read_only() const
{
	return readonly;
}
