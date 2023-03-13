/* *************************************************************

Copyright 2019 Terry Lyons.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */

// libalgebra functionality
#include "libalgebra/alg_types.h"

// boost dependencies for current tests
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp> // exists
#include <boost/filesystem/fstream.hpp>

// std:: dependencies for current tests
#include <iostream>
#include <vector>


// the unit test framework
#include <UnitTest++.h>

// a debugging tool - SHOW(X) outputs variable name X and its content to a stream (e.g. cout) 
#include "SHOW.h"
#include "../common/time_and_details.h"
#include "../common/rng.h"


namespace {
	// the determining template variables
	const unsigned DEPTH = 9;
	const unsigned ALPHABET_SIZE = 2;
	const unsigned STEPS = 1000;

	/// brown_path - a class exposing a fixed Brownian path of dimension "ALPHABET_SIZE" 
	/// the path is evenly sampled at times 0, 1/Steps, 2/Steps, ..., 1.0 within the 
	/// unit interval. The path is presented as a one dimensional stream of doubles [begin(),end())
	struct brown_path
	{
		// state
		const size_t steps;
		const size_t width;
		std::vector<double> path;

		// accessors
			// dimensions
		size_t Steps() const { return steps; }
		size_t Width() const { return width; }

		// the path data
		const double* begin() const { return &(*path.begin()); }
		const double* end() const { return begin() + path.size(); }

		// constructor
		brown_path() :
			steps(STEPS),
			width(ALPHABET_SIZE),
			path((Steps() + 1)* Width(), 0.)
		{
			// set up random number generation
				// seeded generator	
			unsigned int seed = (const unsigned int&)0x6d35f0e5b8f6c603;//std::random_device seed; unsigned int seed = seed();
			mt19937 generator;
			generator.seed(seed);

			// rng
			NORMAL_DIST<double> distribution(0., 1. / sqrt(Steps()));//distribution(mean, std deviation)

			// create random path with dimension width and steps increments, and so a (steps+1) x width C matrix. 
			for (size_t i = 0; i < Steps(); ++i)
				for (size_t j = 0; j < Width(); ++j) {
					double increment = distribution(generator);
					path[(std::size_t(i) + 1) * Width() + j] = path[std::size_t(i) * width + j] + increment;
				}
		}

		// destructor
		~brown_path() {}
	};

	/// brown_path_increments - This class is intended to be derived from in unit tests.
	/// It exposes a fixed Brownian path of dimension "ALPHABET_SIZE", and the path is 
	/// presented via the LIE representation of the (truncated to level "DEPTH") of the 
	/// path increments over the intervals 
	///         [0, 1/Steps), [1/Steps, 2/Steps), ..., [(Steps-1)/Steps, 1.0)
	/// within the unit interval. The path is presented as a stream of size "STEPS" of LIE elements 
	/// in a range [increments.begin(),increments.end())
	///
	/// In addition to the LIE representation of the path increments it provides the basic tensor objects and the maps 
	/// as well as helper functions that can calculate signatures and log signatures of the path over subintervals [b,e)
	///
	struct brown_path_increments : private brown_path
	{
		// types
		typedef alg_types<DEPTH, ALPHABET_SIZE, SPReal>  ALG;
		typedef typename ALG::TENSOR TENSOR;
		typedef typename ALG::LIE LIE;
		typedef typename ALG::MAPS MAPS;
		typedef typename ALG::CBH CBH;
		typedef typename ALG::S S;

		// state
		mutable MAPS maps; // has cached state
		mutable CBH cbh; // has cached state
		using brown_path::Steps;
		using brown_path::Width;
		std::vector<LIE> increments;

		// helper functions
			/// takes flat vector to increment of LIE type
		template<class const_scalar_iterator>
		LIE LieDifference(const_scalar_iterator i) const
		{
			LIE increment;
			for (size_t j = 0; j < Width(); j++)
				// Hall basis elements start at index 1 with zero as a reserved parent index		
				increment += LIE(j + 1, (S) (*(i + j + Width()) - *(i + j)));
			return increment;
		}

		/// brown_path_increments inherits its constructor
		template<typename ...Ts>
		brown_path_increments(Ts... args)
			: brown_path(args...) {
			// convert sampled path to LIE increments
			typename std::vector<double>::const_iterator iend (path.end());
			for (typename std::vector<double>::const_iterator i = path.begin(); (i + Width()) != iend; i += Width()) {
				increments.push_back(LieDifference(i));
			}
			// remove path data
			brown_path::path.clear();
		}

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
			TENSOR signature = signature(begin, end);
			return maps.t2l(log(signature));
		}
	};

	TEST_FIXTURE(brown_path_increments, speed)
	{
		TEST_DETAILS();
		TENSOR sig;
		for (int count = 0; count < 100; ++count)
			signature(increments.begin(), increments.end());
	}

} // namespace
