#pragma once
#include "alg_framework.h"
#include "makebm.h"
#include <vector>

	template <unsigned DEPTH, unsigned ALPHABET_SIZE, unsigned STEPS>
	struct brown_path_increments : alg_framework<DEPTH, ALPHABET_SIZE, DPReal>
	{
	
		typedef typename alg_framework<DEPTH, ALPHABET_SIZE, DPReal>::LIE LIE;
		std::vector<LIE> increments;

		// generates the increments of a brownian path as LIE increments
		brown_path_increments(const unsigned steps = STEPS)
		{
			std::vector<double> path;
			path.reserve(steps);
			makebm(path, steps, ALPHABET_SIZE);
			// convert sampled path to LIE increments
			typename std::vector<double>::const_iterator pend(path.end());
			for (typename std::vector<double>::const_iterator i = path.begin(); (i + ALPHABET_SIZE) != pend; i +=
			        ALPHABET_SIZE) {
				increments.push_back(LieDifference(i, ALPHABET_SIZE));
			}
		}

		typedef typename alg_framework<DEPTH, ALPHABET_SIZE, DPReal>::TENSOR TENSOR;
		typedef typename alg_framework<DEPTH, ALPHABET_SIZE, DPReal>::S S;
		using alg_framework<DEPTH, ALPHABET_SIZE, DPReal>::maps;
		using alg_framework<DEPTH, ALPHABET_SIZE, DPReal>::cbh;
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

	private:
		template<class const_scalar_iterator>
		LIE LieDifference(const_scalar_iterator i, unsigned width) const
		{
			LIE increment;
			for (size_t j = 0; j < width; j++)
				// Hall basis elements start at index 1 with zero as a reserved parent index		
				increment += LIE(j + 1, *(i + j + width) - *(i + j));
			return increment;
		}
	};



