#pragma once
// the libalgebra framework
#include "alg_framework.h"
#include "../common/rng.h"

template <typename alg::LET depth, typename alg::LET width, coefficient_t  S_t = Rational>
struct sigtools: alg_framework <depth, width, S_t>
{
	typedef alg_framework <depth, width, S_t> ALG_FRAMEWORK;
	typedef typename ALG_FRAMEWORK::TENSOR TENSOR;
	typedef typename ALG_FRAMEWORK::LIE LIE;
	typedef typename ALG_FRAMEWORK::S S;
	using ALG_FRAMEWORK::maps;
	using ALG_FRAMEWORK::cbh;

// helper functions

	/// computes a signature from an iterable sequence of lie elements
	template<typename ITERATOR_T>
	TENSOR signature(ITERATOR_T begin, ITERATOR_T end) const
	{
		ptrdiff_t N = end - begin;
		TENSOR signature(S(1));

		for (ptrdiff_t i = 0; i < N; i++)
			signature *= exp(maps.l2t(*(begin + i)));
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

template <typename alg::LET DEPTH, typename alg::LET ALPHABET_SIZE, coefficient_t  S_t = Rational>
struct categorical_path : sigtools<DEPTH, ALPHABET_SIZE, S_t>
{
	typedef sigtools<DEPTH, ALPHABET_SIZE, S_t> SIGTOOLS;
	typedef typename SIGTOOLS::TENSOR TENSOR;
	typedef typename SIGTOOLS::LIE LIE;
	typedef typename SIGTOOLS::S S;

	using SIGTOOLS::signature;
	using SIGTOOLS::logsignature;
	using SIGTOOLS::maps;
	using SIGTOOLS::cbh;

	// state
	const size_t steps;
	const size_t width;
	const size_t depth;
	std::vector<LIE> increments;

	// accessors
			// dimensions
	size_t Steps() const { return steps; }
	size_t Width() const { return width; }
	size_t Depth() const { return depth; }

	// the path data
	const LIE* begin() const { return &(*increments.begin()); }
	const LIE* end() const { return &(*increments.begin()) + increments.size(); }

	// constructor
	categorical_path(size_t STPS = ALPHABET_SIZE) :
		steps(STPS), //number of increments in walk
		width(ALPHABET_SIZE), //constant set in alg_types
		depth(DEPTH), //constant set in alg_types
		increments(Steps(), LIE())
	{
		// set up random number generation
			// seeded generator	
		unsigned int seed = (const unsigned int&)0x6d35f0e5b8f6c603;//std::random_device seed; unsigned int seed = seed();
		mt19937 generator;
		generator.seed(seed);

		// rng
		UNIFORM_INT_DIST<alg::LET> distribution(1, ALPHABET_SIZE);

		// create random categorical increments
		for (size_t i = 0; i != Steps(); ++i)
		//	increments[i] = LIE(distribution(generator), S(1));
			increments[i] = LIE(1+generator() %ALPHABET_SIZE, S(1));
	}
};

// ? OMP
//https://software.intel.com/en-us/node/583439
//page 180 of: http://www.openmp.org/mp-documents/OpenMP4.0.0.pdf		
//#pragma omp declare reduction (mymult:TENSOR:omp_out *= exp(maps.l2t(*(begin + omp_in)))) initializer(omp_priv=TENSOR(S(1)))
//#pragma omp parallel for reduction(mymult:signature)
