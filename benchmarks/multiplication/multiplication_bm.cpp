#include <benchmark/benchmark.h>
#include <iostream>

//#include <libalgebra/libalgebra.h>
//#include <libalgebra/giles_multiplication.h>

template<int W, int D>
static void BM_multiplication(benchmark::State& state)
{
    std::cout << "Hello, Benchmarks!" << std::endl;

    for (auto _ : state) {
////        auto result = tiled_inverse(t);
//        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
//
//    state.SetBytesProcessed(state.iterations()*2 * sizeof(double) * playground::tensor_alg_size(W, D));
//    state.counters["block_letters"] = tiled_inverse.block_letters;
//    state.counters["block_size"] = tiled_inverse.block_size * sizeof(double);
//    state.counters["blob size"] = sizeof(double) * playground::tensor_alg_size(W, D);

}

BENCHMARK_TEMPLATE(BM_multiplication, 1, 1);