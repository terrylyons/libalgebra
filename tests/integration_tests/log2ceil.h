#pragma once

template<typename UNSIGNEDINT>
UNSIGNEDINT log2ceil(UNSIGNEDINT n)
{
	assert((!std::numeric_limits<UNSIGNEDINT>::is_signed) && (std::numeric_limits<UNSIGNEDINT>::digits <= 64));
	assert(n != 0);
	--n;
	n |= n >> 1;
	n |= n >> 2;
	n |= n >> 4;
	n |= n >> 8;
	n |= n >> 16;
#pragma warning (disable : 4293)
	if (std::numeric_limits<UNSIGNEDINT>::digits > 32)
		n |= n >> 32;
#pragma warning (default : 4293)
	++n;
	return n;
}

