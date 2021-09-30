#include "makebm.h"
#include <vector>
#include "../common/rng.h"

void makebm(std::vector<double>& pathi, const size_t steps, const size_t width)
{

	// set up random number generation
	unsigned int seed = (const unsigned int&)0x6d35f0e5b8f6c603;//std::random_device seed; unsigned int seed = seed();
	mt19937 generator;
	generator.seed(seed);
	NORMAL_DIST<double> distribution(0., 1. / sqrt(steps));//distribution(mean, std deviation)

	// create random path with dimension width and steps increments, and so a (steps+1) x width C matrix. 
	std::vector<double> path((steps + 1) * width, 0.);
	for (size_t i = 0; i < steps; ++i)
		for (size_t j = 0; j < width; ++j) {
			double increment = distribution(generator);
			path[(std::size_t(i) + 1) * width + j] = path[std::size_t(i) * width + j] + increment;
		}

	// return it
	pathi.swap(path);
}