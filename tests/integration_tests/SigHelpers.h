#pragma once
#include "TreeBufferHelper.h"
#include <vector>
#include <omp.h>
#include "SHOW.h"
#include "log2ceil.h"
// these helper functions are used widely in the tests
// however, the tests duplicate the code and need to be modified to use these standard versions
// these standard versions have not been tested



/// computes a signature from an iterable sequence of lie elements
template<typename ITERATOR_T, typename FRAMEWORK>
typename FRAMEWORK::TENSOR signature(ITERATOR_T begin, ITERATOR_T end, const FRAMEWORK& context)
{
	// omp friendly
	ptrdiff_t N = end - begin;
	typename FRAMEWORK::TENSOR signature(typename FRAMEWORK::S(1));
	for (ptrdiff_t i = 0; i < N; i++)
		signature *= exp(context.maps.l2t(*(begin + i)));
	return signature;
}

/// computes the logsignature from the signature
template<class ITERATOR_T, typename FRAMEWORK>
typename FRAMEWORK::LIE logsignature(ITERATOR_T begin, ITERATOR_T end, const FRAMEWORK& context)
{
	typename FRAMEWORK::TENSOR sig = signature(begin, end, context);
	return context.maps.t2l(context.log(sig));
}

/// computes the logsignature from the signature using omp
template<class ITERATOR_T, typename FRAMEWORK>
typename FRAMEWORK::LIE o_logsignature(ITERATOR_T begin, ITERATOR_T end, const FRAMEWORK& context)
{
	typename FRAMEWORK::TENSOR sig = signature(begin, end, context);
	return context.maps.t2l(context.log(sig));
}

/// computes a signature from an iterable sequence of lie elements using OMP
template<typename ITERATOR_T, typename FRAMEWORK>
typename FRAMEWORK::TENSOR o_signature(ITERATOR_T begin, ITERATOR_T end, const FRAMEWORK& context)
{
	typedef typename FRAMEWORK::TENSOR TENSOR;
	typedef typename FRAMEWORK::LIE LIE;
	typedef typename FRAMEWORK::S S;
	typename FRAMEWORK::MAPS& maps = context.maps;
	//typename FRAMEWORK::CBH& cbh = context.cbh;
	//#undef  _OPENMP
//#ifndef _OPENMP
#if 0
		 //simple non-parallel form
	TENSOR signature(S(1));
	for (ITERATOR_T i = begin; i != end; i++)
		signature *= exp(maps.l2t(*i));
	return signature;
#else
	ptrdiff_t leaves = log2ceil(size_t(end - begin));
	const LIE* in = &(*begin);
	CTreeBufferHelper t(1, leaves);
	ptrdiff_t e = t.end();
	// pack to power of two with the identity so that helper functions gets order correct
	// storage for the intermediate calculations 
	// only requires one element per thread; 
	// the approach here retains all intermediate values 
	// and is wasteful unless required
	std::vector<TENSOR> sigs(e, TENSOR(S(1)));

#pragma omp parallel for
	for (ptrdiff_t i = 0; i < end - begin; i++)
		sigs[i] = exp(maps.l2t(in[i]));

	for (ptrdiff_t j = t.parent(0); j < e; j = t.parent(j))
	{
		ptrdiff_t jj = t.parent(j);
#pragma omp parallel for
		for (ptrdiff_t i = j; i < jj; i++) {
			// use move semantics and avoid a copy?
			sigs[i] = sigs[t.left(i)] * sigs[t.right(i)];		
		}
	}
	return *(sigs.rbegin());
#endif // _OPENMP
}
