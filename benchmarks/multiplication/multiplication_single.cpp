//
// Created by user on 16/01/23.
//

#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <libalgebra/libalgebra.h>
#include <benchmark/benchmark.h>

#include <thread>

#ifndef LA_TM_BM_WIDTH
#define LA_TM_BM_WIDTH 4
#endif

#ifndef LA_TM_BM_DEPTH
#define LA_TM_BM_DEPTH 13
#endif

#ifndef LA_TM_BM_TILE_LETTERS
#define LA_TM_BM_TILE_LETTERS 3
#endif

#ifndef LA_TM_BM_CACHE_LETTERS
#define LA_TM_BM_CACHE_LETTERS 0
#endif

namespace po = boost::program_options;


template <unsigned W, unsigned D>
using multiplication = alg::tiled_free_tensor_multiplication<W, D, LA_TM_BM_TILE_LETTERS, LA_TM_BM_CACHE_LETTERS>;

using traits = alg::dtl::multiplication_traits<multiplication<LA_TM_BM_WIDTH, LA_TM_BM_DEPTH>>;

using tensor = alg::free_tensor<alg::coefficients::float_field, LA_TM_BM_WIDTH, LA_TM_BM_DEPTH,
      alg::vectors::dense_vector, multiplication>;

void tm_benchmark(int warmup, int iterations) {
    alg::tensor_basis<LA_TM_BM_WIDTH, LA_TM_BM_DEPTH>::KEY kunit;
    tensor lhs, rhs, result;
    lhs.base_vector().resize_to_degree(LA_TM_BM_DEPTH);
    lhs.base_vector().construct_reverse_data(LA_TM_BM_DEPTH - 1);
    lhs[kunit] = 0.3453453f;
    rhs.base_vector().resize_to_degree(LA_TM_BM_DEPTH);
    rhs[kunit] = 0.7324534f;
    result.base_vector().resize_to_degree(LA_TM_BM_DEPTH);

    multiplication<LA_TM_BM_WIDTH, LA_TM_BM_DEPTH> mul;

    for (int i = 0; i < warmup; ++i) {
        benchmark::DoNotOptimize(result.base_vector().as_mut_ptr());
        traits::multiply_and_add(mul, result, lhs, rhs);
        benchmark::ClobberMemory();
    }

    {
//        boost::timer::auto_cpu_timer timer;
        for (int i = 0; i < iterations; ++i) {
            benchmark::DoNotOptimize(result.base_vector().as_mut_ptr());
            result = lhs * rhs;
            benchmark::ClobberMemory();
        }
    }
}


int main(int argc, char** argv) {

    int warmup, iterations;
    po::options_description descr;
    descr.add_options()
            ("warmup", po::value<int>(&warmup)->default_value(0), "warmup cycles")
            ("repetitions", po::value<int>(&iterations)->default_value(5), "timed repetitions");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, descr), vm);
    po::notify(vm);

    std::thread run_bm(tm_benchmark, warmup, iterations);
    run_bm.join();


    return 0;
}
