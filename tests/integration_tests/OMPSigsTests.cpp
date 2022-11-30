
// the libalgebra framework
#include "alg_framework.h"

// the unit test framework
#include <UnitTest++.h>
#include "../common/time_and_details.h"

// omp
#include <vector>
#include "brown_path_increments.h"
#include "SigHelpers.h"
#include <iostream>

// validates the hall set and lie multiplication over it
SUITE(OMP_SIGNATURES)
{
	// DEPTH, ALPHABET SIZE, STEPS
	typedef brown_path_increments<6, 5, 50> SETUP65;
	typedef brown_path_increments<1, 5, 50> SETUP15;
	typedef brown_path_increments<2, 2, 50> SETUP22;
	TEST_FIXTURE(SETUP65, bm65_test_parallel_signature)
	{
		TEST_DETAILS();
		TENSOR sig, sig1, sig2;
		std::cout << "framework bound sequential signature: ";
		{
			timer seqsig_t;
			sig = signature(increments.begin(), increments.end());
		}
		std::cout << "template sequential signature: ";
		{
			timer seqsig_t;
			sig1 = ::signature(increments.begin(), increments.end(), *this);
		}
		std::cout << "parallel signature: ";
		{
			timer seqsig_t;
			sig2 = ::o_signature(increments.begin(), increments.end(), *this);
		} 
		//std::cout << sig << " \n" << sig1 << std::endl;
		//std::cout << sig.NormL1() << " " << (sig - sig12).NormL1() << std::endl;
		CHECK_CLOSE(0., (sig - sig1).NormL1() + (sig2 - sig1).NormL1(), 10e-13);
	}
}
