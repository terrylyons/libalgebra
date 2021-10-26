#pragma once

// the libalgebra framework
#include <libalgebra/libalgebra.h>
#include <libalgebra/alg_types.h>

// simple framework for using libalgebra
template <unsigned DEPTH, unsigned ALPHABET_SIZE, enum coefficient_t scalar_t>
struct alg_framework : alg_types<DEPTH, ALPHABET_SIZE, scalar_t>
{
	// lib algebra required state
	mutable typename alg_types<DEPTH, ALPHABET_SIZE, scalar_t>::MAPS maps;
	mutable typename alg_types<DEPTH, ALPHABET_SIZE, scalar_t>::CBH  cbh;
};

