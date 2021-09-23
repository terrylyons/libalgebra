#pragma once



#if __cplusplus >= 201103L
#include <chrono>
#else
#include <boost/random.hpp>
#endif


#if __cplusplus >= 201103L
namespace chrono = std::chrono;
#else
namespace chrono = boost::chrono;
#endif

struct timer
{
	chrono::time_point<chrono::steady_clock> start;
	chrono::time_point<chrono::steady_clock> stop;

	timer() :start(chrono::steady_clock::now())
	{}

	~timer()
	{
		Adaptive();
	}

	void Adaptive()
	{
		stop = chrono::steady_clock::now();
		chrono::nanoseconds time = stop - start;

		long long ntime;
		if ((ntime = chrono::duration_cast<chrono::nanoseconds>(time).count())< 10000LL)
		std::cout << "Elapsed time : "
			<< ntime
			<< " ns" << std::endl;
		else if ((ntime = chrono::duration_cast<chrono::microseconds>(time).count()) < 10000LL)
			std::cout << "Elapsed time : "
			<< ntime
			<< " micros" << std::endl;
		else if ((ntime = chrono::duration_cast<chrono::milliseconds>(time).count()) < 10000LL)
			std::cout << "Elapsed time : "
			<< ntime
			<< " ms" << std::endl;
		else {
			ntime = chrono::duration_cast<chrono::seconds>(time).count();
			std::cout << "Elapsed time : "
				<< ntime
				<< " s" << std::endl;
		}
	}

	void NanoSecs()
	{
		stop = chrono::steady_clock::now();
		std::cout << "Elapsed time in nanoseconds : "
			<< chrono::duration_cast<chrono::nanoseconds>(stop - start).count()
			<< " ns" << std::endl;
	}

	void MicroSecs()
	{
		stop = chrono::steady_clock::now();
		std::cout << "Elapsed time in microseconds : "
			<< chrono::duration_cast<chrono::microseconds>(stop - start).count()
			<< " �s" << std::endl;
	}

	void MilliSecs()
	{
		stop = chrono::steady_clock::now();
		std::cout << "Elapsed time in milliseconds : "
			<< chrono::duration_cast<chrono::milliseconds>(stop - start).count()
			<< " ms" << std::endl;
	}

	void Secs()
	{
		stop = chrono::steady_clock::now();
		std::cout << "Elapsed time in seconds : "
			<< chrono::duration_cast<chrono::seconds>(stop - start).count()
			<< " sec";
	}
};
#if 0
// a macro to identify and time a test to cout
#ifndef TEST_DETAILS
#ifndef SHO_NO
#define TEST_DETAILS() std::cout << "\nFile : " << m_details.filename << "\nLine : " << m_details.lineNumber << "\nSuite : " << m_details.suiteName << "\nTest : " << m_details.testName << std::endl ; timer t46390__LINE__;
#else
#define TEST_DETAILS()
#endif
#else
#endif
#else
#define TEST_DETAILS()
#endif