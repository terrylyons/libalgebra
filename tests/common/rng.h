//
// Created by sam on 24/02/2021.
//

#ifndef LIBALGEBRAUNITTESTS_RNG_H
#define LIBALGEBRAUNITTESTS_RNG_H

#if __cplusplus >= 201103L
#include <random>
typedef std::mt19937 mt19937;
#define NORMAL_DIST std::normal_distribution
#define UNIFORM_INT_DIST std::uniform_int_distribution
#else
#include <boost/random.hpp>
typedef boost::mt19937 mt19937;
#define NORMAL_DIST boost::normal_distribution
#define UNIFORM_INT_DIST boost::random::uniform_int_distribution
#endif


#endif //LIBALGEBRAUNITTESTS_RNG_H
