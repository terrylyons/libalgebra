/* Copyright -  All Rights Reserved - Terry Lyons 2008 */
#pragma once
#include <stddef.h>  //ptrdiff_t
#include <algorithm> //min
#ifdef _MSC_VER
#include <crtdbg.h>  //_ASSERT
#else
#define _ASSERT(expr) (void(0))
#endif


class CTreeBufferHelper
{
	// the number of trees in the initial forest
	ptrdiff_t iNoTrees;
	// the number of leaves in the initial forest
	ptrdiff_t iInitialNoLeaves;

public:
	// Constructor
	CTreeBufferHelper(ptrdiff_t SmallestReducibleSetSize, ptrdiff_t NoPointsToBeprocessed)
		: iNoTrees(SmallestReducibleSetSize),
		iInitialNoLeaves(NoPointsToBeprocessed)
	{
		_ASSERT(iInitialNoLeaves >= iNoTrees && iNoTrees > 0);
	}

	bool isleaf(const ptrdiff_t& node) const
	{
		return (node < iInitialNoLeaves && node >= 0);
	}

	ptrdiff_t end() const
	{
		return 2 * iInitialNoLeaves - iNoTrees;
	}

	bool isnode(const ptrdiff_t& node) const
	{
		return node >= 0 && node < end();
	}

	// in the reduction all dependencies of [node,parent) are const in [0,node) 
	// and can be computed in parallel
	ptrdiff_t parent(const ptrdiff_t& node) const
	{
		_ASSERT(isnode(node));
		return std::min(iInitialNoLeaves + (node / 2), end());
	}

	bool isroot(const ptrdiff_t& node) const
	{
		_ASSERT(isnode(node));
		return parent(node) == end();
	}

	ptrdiff_t left(const ptrdiff_t& node) const
	{
		_ASSERT(isnode(node));
		// returns negative if leaf
		return (node - iInitialNoLeaves) * 2;
	}

	ptrdiff_t right(const ptrdiff_t& node) const
	{
		return (left(node) < 0) ? -1 : left(node) + 1;
	}

};
