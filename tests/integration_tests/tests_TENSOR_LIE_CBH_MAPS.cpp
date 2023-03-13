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
#include "../common/rng.h"

// to allow redefinitions in other test modules
namespace {
	// the tensor shape variables
	const unsigned DEPTH = 5;
	const unsigned ALPHABET_SIZE = 5;
	const unsigned DEPTH2 = 4;
	const unsigned ALPHABET_SIZE2 = 3;
	const unsigned STEPS = 30;

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

	struct categorical_path : uni_env2
	{
		// state
		const size_t steps;
		const size_t width;
		std::vector<LIE> increments;

		// accessors
				// dimensions
		size_t Steps() const { return steps; }
		size_t Width() const { return width; }

		// the path data
		const LIE* begin() const { return &(*increments.begin()); }
		const LIE* end() const { return begin() + increments.size(); }

		// constructor
		categorical_path() :
			steps(STEPS),
			width(ALPHABET_SIZE),
			increments((Steps()), LIE())
		{
			// set up random number generation
				// seeded generator	
			unsigned int seed = (const unsigned int&)0x6d35f0e5b8f6c603;//std::random_device seed; unsigned int seed = seed();
			mt19937 generator;
			generator.seed(seed);

			// rng
			UNIFORM_INT_DIST<LET> distribution(1, ALPHABET_SIZE);

			// create random categorical increments
			for (size_t i = 0; i < Steps(); ++i)
				increments[i] = LIE(distribution(generator), S(1));
		}
		// helper functions
				/// computes a signature from an iterable sequence of lie elements
		template<class ITERATOR_T>
		TENSOR signature(ITERATOR_T begin, ITERATOR_T end) const
		{
			TENSOR signature(S(1));
			for (ITERATOR_T i = begin; i != end; i++)
				signature *= exp(maps.l2t(*i));
			return signature;
		}

		/// computes the logsignature from the signature
		template<class ITERATOR_T>
		LIE logsignature(ITERATOR_T begin, ITERATOR_T end) const
		{
			TENSOR sig = signature(begin, end);
			return maps.t2l(log(sig));
		}
	};

	TEST_FIXTURE(categorical_path, long_multiplication)
	{
		TEST_DETAILS();
		typename std::vector<LIE>::const_iterator begin = increments.begin();
		typename std::vector<LIE>::const_iterator end = increments.end();
		TENSOR sig = signature(begin, end);
		for (typename std::vector<LIE>::const_iterator i = begin; i != end; i++) {
			TENSOR err = sig - signature(begin, i) * signature(i, end);
			CHECK_EQUAL(TENSOR(), err);
		}
	}

	TEST_FIXTURE(categorical_path, bug)
	{
		TEST_DETAILS();
		//4: { 1(4) }
		//446 : { 1([3, [2, [2, [2, 3]]]]) }

		LIE l1(1), l2(2), l3(3), l4(4);
		LIE k1(l4);
		LIE k2(l3 * (l2 * (l2 * (l2 * (l3)))));
		LIE ans = k1 * k2;
	}

	TEST_FIXTURE(categorical_path, units)
	{
		TEST_DETAILS();
		LET let(2);
		S scalar(2);

		LIE lzero;
		//LIE lscalar(scalar); // not defined
		LIE lnotscalar(let);

		TENSOR tzero;
		TENSOR tone(S(1));

		TENSOR tscalar(scalar);
		TENSOR tnotscalar(let, S(1));

		// note LIE(S(4)) should be undefined behavior
		// LIE(4) is not but treats the argument as a letter
		typename std::vector<LIE>::const_iterator begin = increments.begin();
        typename std::vector<LIE>::const_iterator end = increments.end();
		TENSOR sig = signature(begin, end);
		LIE logsig = logsignature(begin, end);

		CHECK_EQUAL(logsig, logsig + lzero);
		CHECK_EQUAL(logsig, lzero + logsig);
		CHECK_EQUAL(lzero, lzero + lzero);

		CHECK_EQUAL(lzero, logsig * lzero);
		CHECK_EQUAL(lzero, lzero * logsig);
		CHECK_EQUAL(lzero, lzero * lzero);

		CHECK_EQUAL(sig, sig + tzero);
		CHECK_EQUAL(sig, tzero + sig);
		CHECK_EQUAL(tzero, tzero + tzero);

		CHECK_EQUAL(tzero, sig * tzero);
		CHECK_EQUAL(tzero, tzero * sig);
		CHECK_EQUAL(tzero, tzero * tzero);

		CHECK_EQUAL(sig, sig * tone);
		CHECK_EQUAL(sig, tone * sig);
		CHECK_EQUAL(tone, tone * tone);

		CHECK_EQUAL(sig * scalar, sig * tscalar);
		CHECK_EQUAL(sig * scalar, tscalar * sig);
		CHECK_EQUAL(tone * scalar * scalar, tscalar * tscalar);

		// danger
		CHECK(LIE(2) == LIE(2, 1));
		CHECK(LIE(0) != LIE());
		// 0 is not a letter in the LIE basis
		CHECK(tnotscalar != tscalar);
		// danger!! since there are conversions from int types to scalar
	}

	TEST_FIXTURE(categorical_path, CBH)
	{
		TEST_DETAILS();
        typename std::vector<LIE>::const_iterator begin = increments.begin();
        typename std::vector<LIE>::const_iterator end = increments.end();
		TENSOR sig = signature(begin, end);
		LIE logsig = logsignature(begin, end);
		CHECK(((exp(maps.l2t(logsig)) - sig) == TENSOR()));
		CHECK(maps.t2l(log(sig)) == logsig);
		CHECK(exp(log(sig)) == sig);
	}

	TEST_FIXTURE(categorical_path, LIE_PRODUCT)
	{
		TEST_DETAILS();
        typename std::vector<LIE>::const_iterator begin = increments.begin();
        typename std::vector<LIE>::const_iterator end = increments.end();
		size_t midpoint = (end - begin) / 2;
		LIE l1 = logsignature(begin, begin + midpoint);
		LIE l2 = logsignature(begin + midpoint, end);
		TENSOR t1 = maps.l2t(l1);
        TENSOR t2 = maps.l2t(l2);
		TENSOR t = t1 * t2 - t2 * t1;
		CHECK((maps.l2t(l1 * l2) - t) == TENSOR());
	}

	TEST_FIXTURE(categorical_path, AntipodeMap)
	{
		TEST_DETAILS();
        typename std::vector<LIE>::const_iterator begin = increments.begin();
        typename std::vector<LIE>::const_iterator end = increments.end();
		LIE l1 = logsignature(begin, end);
		TENSOR t1 = maps.l2t(l1);

		CHECK(inverse(exp(t1)) == reflect(exp(t1)));
		CHECK(inverse(exp(t1))* exp(t1) == TENSOR(S(1)));
	}

	TEST_FIXTURE(uni_env, the_basis)
	{
		TEST_DETAILS();
		// show the bases:
		{
			std::string expected(" 1 2 [1,2] [1,[1,2]] [2,[1,2]] [1,[1,[1,2]]] [2,[1,[1,2]]] [2,[2,[1,2]]] [1,[1,[1,[1,2]]]] [2,[1,[1,[1,2]]]] [2,[2,[1,[1,2]]]] [2,[2,[2,[1,2]]]] [[1,2],[1,[1,2]]] [[1,2],[2,[1,2]]]");
			std::string ans;
			for (typename LIE::BASIS::KEY k = LIE::basis.begin(); k != LIE::basis.end();
				k = LIE::basis.nextkey(k))
				ans += std::string(" ") + LIE::basis.key2string(k);
			//CHECK_EQUAL(expected, ans);
		}
		{
			std::string expected(" () (1) (2) (1,1) (1,2) (2,1) (2,2) (1,1,1) (1,1,2) (1,2,1) (1,2,2) (2,1,1) (2,1,2) (2,2,1) (2,2,2) (1,1,1,1) (1,1,1,2) (1,1,2,1) (1,1,2,2) (1,2,1,1) (1,2,1,2) (1,2,2,1) (1,2,2,2) (2,1,1,1) (2,1,1,2) (2,1,2,1) (2,1,2,2) (2,2,1,1) (2,2,1,2) (2,2,2,1) (2,2,2,2) (1,1,1,1,1) (1,1,1,1,2) (1,1,1,2,1) (1,1,1,2,2) (1,1,2,1,1) (1,1,2,1,2) (1,1,2,2,1) (1,1,2,2,2) (1,2,1,1,1) (1,2,1,1,2) (1,2,1,2,1) (1,2,1,2,2) (1,2,2,1,1) (1,2,2,1,2) (1,2,2,2,1) (1,2,2,2,2) (2,1,1,1,1) (2,1,1,1,2) (2,1,1,2,1) (2,1,1,2,2) (2,1,2,1,1) (2,1,2,1,2) (2,1,2,2,1) (2,1,2,2,2) (2,2,1,1,1) (2,2,1,1,2) (2,2,1,2,1) (2,2,1,2,2) (2,2,2,1,1) (2,2,2,1,2) (2,2,2,2,1) (2,2,2,2,2)");
			std::string ans;
			for (typename TENSOR::BASIS::KEY k = TENSOR::basis.begin();
				k != TENSOR::basis.end(); k = TENSOR::basis.nextkey(k))
				ans += std::string(" (") + TENSOR::basis.key2string(k) +
				std::string(")");
			//CHECK_EQUAL(expected, ans);

		}
	}
	
	TEST_FIXTURE(uni_env, constructors)
	{
		TEST_DETAILS();
		// construct  
		   // a scalar
		S s(3);
		// letters
		LET l(1), r(2);
		// a lie element
		LIE lie, liel(l,S(1)), liels(l, s);
		LIE lielr = liels * LIE(r, S(1));
		// a tensor (constructible on a scalar)
		TENSOR ten, tens(s), tenls(l, s), tens3 = tenls * TENSOR(r, S(1)) * tenls;
		std::stringstream ans_stream;
		ans_stream << l << " " << lie << " " << liel << " " << liels << " " << lielr << " " << ten << " " << tens << " " << tenls << " " << tenls << " " << tens3;
		std::string expected("1 { } { 1(1) } { 3(1) } { 3([1,2]) } { } { 3() } { 3(1) } { 3(1) } { 9(1,2,1) }");
		std::string actual = ans_stream.str();
		CHECK_EQUAL(expected, actual);
	}

	struct uni_env3 : uni_env2
	{
		template <class LIE_T>
		LIE_T adj(const LIE_T& inner) {
			return LIE_T(inner);
		}

		template<class LIE_T, typename ...Ts>
		LIE_T adj(const LIE_T& outer, Ts... args) {
			return outer * adj(args...);
		}
	};

	TEST_FIXTURE(uni_env3, LibAlgebraTestLieandCBH)
	{
		TEST_DETAILS();
		using std::cout;
		using std::endl;
		using std::vector;

		LIE l1(1, 1), l2(2, 3);
		cout << "xxx" << (((l2 *= (S)4) *= (l1 * l2 + l1)), l2 = l2 * l2, (l2 += l1)) << endl;
		cout << maps.l2t(l2) << endl;

		//l2 = LIE(2, 3);
		//TENSOR t1(maps.l2t(l1)), t2(maps.l2t(l2)), t(t1 + t2 + t1 * t2 - t2 * t1);
		//cout << "t: " << t << endl;
		//cout << maps.t2l(t) << " L1 norm is " << (maps.t2l(t)).NormL1() << endl;
		//cout << maps.t2l(t) << " L1 norm of the second order component is " << (maps.t2l(t)).NormL1(2) << endl;

		//cout << "l2t(t2l) test: " << (t == maps.l2t(maps.t2l((t)))) << endl;
		//cout << "log(exp) test: " << (t - log(exp(t))) << " " << (t == log(exp(t))) << endl;
		//cout << "inverse test1: " << (inverse(exp(-t)) - exp(t)) << " " << (inverse(exp(-t)) == exp(t)) << endl;
		//{
		//	int f(2);
		//	cout << "inverse test2: " << inverse(exp(t) * TENSOR(S(f))) * (exp(t) * TENSOR(S(f))) - TENSOR(S(1))
		//		<< " " << ((inverse(exp(t) * TENSOR(S(2))) * (exp(t) * TENSOR(S(2)))) == TENSOR(S(1)))
		//		<< endl;
		//}
		//cout << "reflect test: " << (inverse(exp(t)) - reflect(exp(t))) << " " << (inverse(exp(t)) == reflect(exp(t))) << endl;
		//vector<LET> letters;
		//letters.push_back(1);
		//letters.push_back(2);

		//cout << "Basic CBH test : " << cbh.basic(letters) << endl;
		//LIE a(1, 1);
		//LIE b(2, 1);
		//vector<const LIE*> v;
		//v.push_back(&a);
		//v.push_back(&b);
		//cout << "Full CBH test : ";
		//cout << (cbh.full(v) == cbh.basic(letters)) << endl;

		//v.clear();
		//v.push_back(&l1);
		//v.push_back(&l2);
		//cout << "Replace in CBH test: ";
		//cout << (replace(cbh.basic(letters), letters, v) == cbh.full(v)) << endl;

		//{
		//	// Compute the coordinate-wise min of two numerical tensors - tricky since the zero values count
		//	TENSOR tt(t);
		//	tt += TENSOR(2, 4);
		//	TENSOR T(t), TT(tt);
		//	cout << "Compute the coordinate-wise min of two tensors: " << endl;
		//	cout << t << endl;
		//	cout << tt << endl;
		//	cout << (t & tt) << endl;
		//	cout << (tt & t) << endl;
		//	cout << (T &= tt) << endl;
		//	cout << (TT &= t) << endl;

		//	// compute the linear transformation of tensors induced by a map of letters
		//	auto h = [](LET in)->LET {return in + LET(1); };
		//	MAPS2::t2t H(h);
		//	tt *= tt += TENSOR(TENSOR::SCALAR(7));
		//	cout << tt << "\n" << H(tt) << endl;
		//}
	}

}// end of anonymous namespace
