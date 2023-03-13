/* *************************************************************

Copyright 2010-2019 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurkó and Arend Janssen.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */
// THIS TEST SUITE IS DERIVED FROM LIBALGEBRA-DEMO.cpp
// BUT EXPLOITS THE MORE UNIFIED SETUP INTERFACE NOW AVAILABLE
//
// BREAKING CHANGE in CBH which now consumes vector<const LIE*> not 
// vector<LIE*> and there is no natural conversion. Change your code.

// libalgebra functionality
#include "libalgebra/alg_types.h"

// std:: dependencies for current tests
#include <iostream>
#include <vector>

// the unit test framework
#include <UnitTest++.h>

// a debugging tool - SHOW(X) outputs variable name X and its content to a stream (e.g. cout) 
#include "SHOW.h"
#include "../common/time_and_details.h"

const int DEPTH = 24;
const int WIDTH = 2;
typedef alg_types<DEPTH, WIDTH, Rational> framework;

//TEST(powersoftwo)
//{
//
//}

TEST_FIXTURE(framework, big_num)
{
	TEST_DETAILS();
	MAPS maps;
	LIE arg = LIE(1, S(1));// +LIE(2, -2);
	TENSOR expo = exp(maps.l2t(arg));
	std::cout << "The lie basis with depth "<<DEPTH<<" and width "<<WIDTH << " is: " << LIE::basis.size() << "\n";
	std::cout << "The tensor basis with depth " << DEPTH << " and width " << WIDTH << " is: " << TENSOR::basis.size() << "\n";
	size_t tensor_dim = TENSOR::basis.size() + 1;
	size_t pwr2 = alg::ConstPower< 2ULL, DEPTH + 1>::ans;
	CHECK_EQUAL(pwr2, tensor_dim );
}
