/* *************************************************************

Copyright 2010-2019 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurko and Arend Janssen.

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

namespace {

	// the tensor shape variables
	const unsigned DEPTH = 5;
	const unsigned ALPHABET_SIZE = 2;
	const unsigned DEPTH2 = 4;
	const unsigned ALPHABET_SIZE2 = 3;

	// the number field
		// IEEE float - SPReal,
		// IEEE double - DPReal, 
		// GMP/MPIR rational arithmetic - Rational     
	const coefficient_t SCALAR_T = Rational;

	/// uni_env - a class defining a single algebraic environment
	struct uni_env
	{
		// types
		typedef alg_types<DEPTH, ALPHABET_SIZE, SCALAR_T>  ALG;
		typedef typename ALG::TENSOR TENSOR;
		typedef typename ALG::LIE LIE;
		typedef typename ALG::MAPS MAPS;
		typedef typename ALG::CBH CBH;
		typedef typename ALG::S S;
		typedef typename ALG::LET LET;

		// state
		mutable MAPS maps; // has cached state
		mutable CBH cbh; // has cached state
	};

	/// uni_env2 - a class defining a second algebraic environment side by side with the first
	struct uni_env2 : uni_env
	{
		// types
		typedef alg_types<DEPTH2, ALPHABET_SIZE2, SCALAR_T>  ALG2;
		typedef typename ALG2::TENSOR TENSOR2;
		typedef typename ALG2::LIE LIE2;
		typedef typename ALG2::MAPS MAPS2;
		typedef typename ALG2::CBH CBH2;
		typedef typename ALG2::S S2;
		typedef typename ALG2::LET LET2;

		// state
		mutable MAPS2 maps2; // has cached state
		mutable CBH2 cbh2; // has cached state
	};

	TEST_FIXTURE(uni_env, LibAlgebraTestPolyLie)
	{
		TEST_DETAILS();
		//'checks' the keyofletter function via this constructor
		ALG::POLYLIE poly11(LET(2), S(1));
		std::cout << poly11 << std::endl;
		//uses the explicit constructor and then uses the Lie bracket
		ALG::POLYLIE poly1(1, 1, 2);
		ALG::POLYLIE poly2(1, 2, 1);
		ALG::POLYLIE poly4(2, 1, 1);
		std::cout << poly4 << std::endl;
		poly4 *= S(3);
		std::cout << poly4 << std::endl;
		poly4 += poly1;
		ALG::POLYLIE poly3;
		poly3 = poly4 * poly2;
		std::cout << poly4 << " * " << poly2 << std::endl << poly3;
	}

	TEST_FIXTURE(uni_env2, LibAlgebraTestCommutativePoly)
	{
		TEST_DETAILS();

		typedef ALG::MULTIPOLY1 MULTIPOLY1;// one dimensional polynomials
		MULTIPOLY1 p0(S(7)), p1(1, 1), p2(2, 3), p3, x(0), y(2), z;
		p3 = p0 + p2 * p1 * p2 + p1 * p2;
		std::cout << "  p0 = " << p0 << "  p1(1,1) = " << p1 << "  p2(2,3) = " << p2 << "  p0 + p2*p1*p2+p1*p2 = " << p3 <<
			" and its degree is " << p3.degree() << std::endl;
		x = p3 * p3;
		std::cout << "  p3 * p3 = " << x << "  and its degree is " << x.degree() << std::endl;
		y = exp(p1);
		z = MULTIPOLY1::diff(p3 * p3, 2);
		std::map<LET, S> coordinate;
		coordinate[1] = S(100);
		std::cout << "  the exp of p1 = " << y << " and its value at coord1 = 100 is " << y.eval(coordinate) << std::endl;
		std::cout << "p3 * p3 derivative in the 2 co-ordinate is " << z << std::endl;
	}

    struct test_letter_transform
    {
	    template <typename LET>
        LET operator()(LET in) const { return in + LET(1); }
    };

	TEST_FIXTURE(uni_env2, LibAlgebraTestLieandCBH)
	{
		TEST_DETAILS();
		using std::cout;
		using std::endl;
		using std::vector;

		LIE l1(1, S(1));
		LIE l2(2, S(3));

		cout << l1 << " and " << l2 << endl;
		l2 *= (S)4;
		cout << l1 << " and " << l2 << endl;
		l2 *= l1 * l2 + l1;
		l2 = l2 * l2;
		l2 += l1;
		cout << l1 << " and " << l2 << endl;
		cout << maps.l2t(l2) << endl;

		l2 = LIE(2, S(3));
		TENSOR t1(maps.l2t(l1)), t2(maps.l2t(l2)), t(t1 + t2 + t1 * t2 - t2 * t1);
		cout << "t: " << t << endl;
		cout << maps.t2l(t) << " L1 norm is " << (maps.t2l(t)).NormL1() << endl;
		cout << maps.t2l(t) << " L1 norm of the second order component is " << (maps.t2l(t)).NormL1(2) << endl;

		cout << "l2t(t2l) test: " << (t == maps.l2t(maps.t2l((t)))) << endl;
		cout << "log(exp) test: " << (t - log(exp(t))) << " " << (t == log(exp(t))) << endl;
		cout << "inverse test1: " << (inverse(exp(-t)) - exp(t)) << " " << (inverse(exp(-t)) == exp(t)) << endl;
		{
			int f(2);
			cout << "inverse test2: " << inverse(exp(t) * TENSOR(S(f))) * (exp(t) * TENSOR(S(f))) - TENSOR(S(1))
				<< " " << ((inverse(exp(t) * TENSOR(S(2))) * (exp(t) * TENSOR(S(2)))) == TENSOR(S(1)))
				<< endl;
		}
		cout << "reflect test: " << (inverse(exp(t)) - reflect(exp(t))) << " " << (inverse(exp(t)) == reflect(exp(t))) << endl;
		vector<LET> letters;
		letters.push_back(1);
		letters.push_back(2);

		cout << "Basic CBH test : " << cbh.basic(letters) << endl;
		LIE a(1, 1);
		LIE b(2, 1);
		vector<const LIE*> v;
		v.push_back(&a);
		v.push_back(&b);
		cout << "Full CBH test : ";
		cout << (cbh.full(v) == cbh.basic(letters)) << endl;

		v.clear();
		v.push_back(&l1);
		v.push_back(&l2);
		cout << "Replace in CBH test: ";
		cout << (replace(cbh.basic(letters), letters, v) == cbh.full(v)) << endl;

		{
			// Compute the coordinate-wise min of two numerical tensors - tricky since the zero values count
			TENSOR tt(t);
			tt += TENSOR(2, 4);
			TENSOR T(t), TT(tt);
			cout << "Compute the coordinate-wise min of two tensors: " << endl;
			cout << t << endl;
			cout << tt << endl;
			cout << (t & tt) << endl;
			cout << (tt & t) << endl;
			cout << (T &= tt) << endl;
			cout << (TT &= t) << endl;

			// compute the linear transformation of tensors induced by a map of letters

/*
			MAPS2::t2t H(test_letter_transform);
			tt *= tt += TENSOR(TENSOR::SCALAR(7));
			cout << tt << "\n" << H(tt) << endl;
			*/
		}
	}

}//namespace
 //EOF.
