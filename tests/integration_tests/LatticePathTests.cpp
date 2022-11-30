// FastSigsExp1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include "categorical_path.h"

// the unit test framework
#include <UnitTest++.h>

// a timer class
#include "../common/time_and_details.h"
#include <iostream>
#include <algorithm>

// the test cases
typedef categorical_path<2, 3> CPD2W3;
typedef categorical_path<3, 3> CPD3W3;
typedef categorical_path<4, 4> CPD4W4;
typedef categorical_path<5, 5> CPD5W5;
typedef categorical_path<6, 6> CPD6W6;
typedef categorical_path<7, 7> CPD7W7;
typedef categorical_path<8, 8> CPD8W8;
typedef categorical_path<9, 9> CPD9W9;


template<typename LOGS, typename SIG, typename FRAMEWORK>
void report_outcomes(const LOGS& logs, const SIG& sig, const FRAMEWORK & context)
{
	std::cout 
		<< "\nTesting integer multiplication and sparseness:\n"
		<< "\t ALPHABET_SIZE = " << FRAMEWORK::ALPHABET_SIZE
		<< "\t DEPTH = " << FRAMEWORK::DEPTH
		<< "\t INTERVALS = " << FRAMEWORK::ALPHABET_SIZE
		<< "\nsupport for the logsignature: " << logs.size() << "/" << FRAMEWORK::LIE::basis.size()
		<< "\nand for the signature:        " << sig.size() << "/" << FRAMEWORK::TENSOR::basis.size()
		<< "\n\n";
}

SUITE(lattice_paths)
{

	TEST_FIXTURE(CPD3W3, short_lattice_path_high_dimension)
	{
		categorical_path p;
		categorical_path::TENSOR sig;
		categorical_path::LIE logs;
		{
			TEST_DETAILS();
			sig = p.signature(p.begin(), p.end());
			logs = p.logsignature(p.begin(), p.end());
			LIE u = p.maps.t2l(log(sig));
		}
		//categorical_path::LIE logsig = p.maps.t2l(log(sig));
		report_outcomes(logs, sig, *this);
		CHECK_EQUAL(5, logs.size());
		CHECK_EQUAL(10, sig.size());
		CHECK_EQUAL(14, categorical_path::LIE::basis.size());
		CHECK_EQUAL(40, categorical_path::TENSOR::basis.size());
	}

	TEST_FIXTURE(CPD4W4, short_lattice_path_high_dimension)
	{
		categorical_path p;
		categorical_path::TENSOR sig;
		categorical_path::LIE logs;
		{
			TEST_DETAILS();
			sig = p.signature(p.begin(), p.end());
			logs = p.logsignature(p.begin(), p.end());
			LIE u = p.maps.t2l(log(sig));
		}
		//categorical_path::LIE logsig = p.maps.t2l(log(sig));
		report_outcomes(logs, sig, *this);
		CHECK_EQUAL(4, logs.size());
		CHECK_EQUAL(25, sig.size());
		CHECK_EQUAL(90, categorical_path::LIE::basis.size());
		CHECK_EQUAL(341, categorical_path::TENSOR::basis.size());
	}

	TEST_FIXTURE(CPD5W5, short_lattice_path_high_dimension)
	{
		categorical_path p;
		categorical_path::TENSOR sig;
		categorical_path::LIE logs;
		{
			TEST_DETAILS();
			sig = p.signature(p.begin(), p.end());
			logs = p.logsignature(p.begin(), p.end());
			LIE u = p.maps.t2l(log(sig));
		}
		//categorical_path::LIE logsig = p.maps.t2l(log(sig));
		report_outcomes(logs, sig, *this);
		CHECK_EQUAL(268, logs.size());
		CHECK_EQUAL(182, sig.size());
		CHECK_EQUAL(829, categorical_path::LIE::basis.size());
		CHECK_EQUAL(3906, categorical_path::TENSOR::basis.size());
	}

	TEST_FIXTURE(CPD6W6, short_lattice_path_high_dimension)
	{
		categorical_path p;
		categorical_path::TENSOR sig;
		categorical_path::LIE logs;
		{
			TEST_DETAILS();
			sig = p.signature(p.begin(), p.end());
			logs = p.logsignature(p.begin(), p.end());
			LIE u = p.maps.t2l(log(sig));
		}
		//categorical_path::LIE logsig = p.maps.t2l(log(sig));
		report_outcomes(logs, sig, *this);
		CHECK_EQUAL(156, logs.size());
		CHECK_EQUAL(189, sig.size());
		CHECK_EQUAL(9695, categorical_path::LIE::basis.size());
		CHECK_EQUAL(55987, categorical_path::TENSOR::basis.size());
	}


	TEST_FIXTURE(CPD7W7, short_lattice_path_high_dimension)
	{
		categorical_path p;
		categorical_path::TENSOR sig;
		categorical_path::LIE logs;
		{
			TEST_DETAILS();
			sig = p.signature(p.begin(), p.end());
			logs = p.logsignature(p.begin(), p.end());
			LIE u = p.maps.t2l(log(sig));
		}
		//categorical_path::LIE logsig = p.maps.t2l(log(sig));
		report_outcomes(logs, sig, *this);
		CHECK_EQUAL(13521, logs.size());
		CHECK_EQUAL(2942, sig.size());
		CHECK_EQUAL(141280, categorical_path::LIE::basis.size());
		CHECK_EQUAL(960800, categorical_path::TENSOR::basis.size());
	}

	TEST_FIXTURE(CPD8W8, short_lattice_path_high_dimension)
	{
		categorical_path p;
		categorical_path::TENSOR sig;
		categorical_path::LIE logs;
		{
			TEST_DETAILS();
			sig = p.signature(p.begin(), p.end());
			logs = p.logsignature(p.begin(), p.end());
			LIE u = p.maps.t2l(log(sig));
		}
		//categorical_path::LIE logsig = p.maps.t2l(log(sig));
		report_outcomes(logs, sig, *this);
		CHECK_EQUAL(11245, logs.size());
		CHECK_EQUAL(1965, sig.size());
		CHECK_EQUAL(2447592, categorical_path::LIE::basis.size());
		CHECK_EQUAL(19173961, categorical_path::TENSOR::basis.size());
	}

}
