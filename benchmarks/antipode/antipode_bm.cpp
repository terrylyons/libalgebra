#include <benchmark/benchmark.h>
#include <iostream>

#include <libalgebra/libalgebra.h>
#include <libalgebra/tensor.h>
#include <libalgebra/dense_vector.h>
#include "../../tests/common/random_vector_generator.h"
#include "../../tests/common/rng.h"

using namespace alg;
using alg::LET;
typedef alg::coefficients::coefficient_field<float> float_field;
typedef alg::coefficients::coefficient_field<alg::coefficients::rational> rational_field;
using float_dist = std::uniform_real_distribution<float_field::S>;

template<DEG Width, DEG Depth>
static void BM_memcpy(benchmark::State& state) {

    //    TODO: replace random vector generator for dense vectors
    //    auto* ptr = lhs.base_vector().as_mut_ptr();
    //    std::generate(ptr, ptr + lhs.base_vector().dimension(), [this, distribution]() { return distribution(rng);

    typedef alg::free_tensor<float_field, Width, Depth, alg::vectors::dense_vector> TENSOR;

    typedef typename TENSOR::KEY KEY;

    TENSOR input_tensor;

    TENSOR result;

    result.base_vector().resize_to_degree(Depth);
    input_tensor.base_vector().resize_to_degree(Depth);

    void* dest = result.base_vector().as_mut_ptr();
    const void* src = input_tensor.base_vector().as_ptr();

    for (auto _: state) {
        benchmark::DoNotOptimize(result.base_vector().as_mut_ptr());

        std::memcpy(dest, src, sizeof(float)*input_tensor.base_vector().dimension());

        benchmark::ClobberMemory();
    }

    state.SetBytesProcessed(2*sizeof(float)*state.iterations()*input_tensor.base_vector().dimension());

}

template<DEG Width, DEG Depth, DEG TileLetters>
static void BM_antipode(benchmark::State& state) {

    //    TODO: replace random vector generator for dense vectors
    //    auto* ptr = lhs.base_vector().as_mut_ptr();
    //    std::generate(ptr, ptr + lhs.base_vector().dimension(), [this, distribution]() { return distribution(rng);

    typedef alg::free_tensor<float_field, Width, Depth, alg::vectors::dense_vector> TENSOR;
//    typedef alg::lie<float_field, Width, Depth, alg::vectors::dense_vector> LIE;
//
//    using rvg_t = la_testing::random_vector_generator<TENSOR, float_dist>;
//    using rvg_l = la_testing::random_vector_generator<LIE, float_dist>;
//
//    const TENSOR tunit;
//    const TENSOR tzero;
//
//    std::mt19937 rngt;
//    rvg_t rvgt;
//
//    std::mt19937 rngl;
//    rvg_l rvgl;

    typedef typename TENSOR::KEY KEY;

//    TENSOR input_tensor = rvgt(rngt);
    TENSOR input_tensor;
    input_tensor.base_vector().resize_to_degree(Depth);

    TENSOR result;

    for (auto _: state) {
        benchmark::DoNotOptimize(result.base_vector().as_mut_ptr());
//        result = antipode(input_tensor);
        dtl::tiled_inverse_operator<Width, Depth, float_field , dtl::default_signer, TileLetters>
                ::apply(input_tensor.base_vector(), result.base_vector());
        benchmark::ClobberMemory();
    }

    state.SetBytesProcessed(2*sizeof(float)*state.iterations()*input_tensor.base_vector().dimension());

}

BENCHMARK_TEMPLATE(BM_memcpy, 4, 4);
BENCHMARK_TEMPLATE(BM_memcpy, 4, 5);
BENCHMARK_TEMPLATE(BM_memcpy, 4, 6);
BENCHMARK_TEMPLATE(BM_memcpy, 4, 7);
BENCHMARK_TEMPLATE(BM_memcpy, 4, 8);
BENCHMARK_TEMPLATE(BM_memcpy, 4, 9);
BENCHMARK_TEMPLATE(BM_memcpy, 4, 10);
BENCHMARK_TEMPLATE(BM_memcpy, 4, 11);
BENCHMARK_TEMPLATE(BM_memcpy, 4, 12);
BENCHMARK_TEMPLATE(BM_memcpy, 4, 13);

BENCHMARK_TEMPLATE(BM_antipode, 4, 4, 1);
BENCHMARK_TEMPLATE(BM_antipode, 4, 4, 2);
BENCHMARK_TEMPLATE(BM_antipode, 4, 4, 3);

BENCHMARK_TEMPLATE(BM_antipode, 4, 5, 1);
BENCHMARK_TEMPLATE(BM_antipode, 4, 5, 2);
BENCHMARK_TEMPLATE(BM_antipode, 4, 5, 3);

BENCHMARK_TEMPLATE(BM_antipode, 4, 6, 1);
BENCHMARK_TEMPLATE(BM_antipode, 4, 6, 2);
BENCHMARK_TEMPLATE(BM_antipode, 4, 6, 3);

BENCHMARK_TEMPLATE(BM_antipode, 4, 7, 1);
BENCHMARK_TEMPLATE(BM_antipode, 4, 7, 2);
BENCHMARK_TEMPLATE(BM_antipode, 4, 7, 3);

BENCHMARK_TEMPLATE(BM_antipode, 4, 8, 1);
BENCHMARK_TEMPLATE(BM_antipode, 4, 8, 2);
BENCHMARK_TEMPLATE(BM_antipode, 4, 8, 3);

BENCHMARK_TEMPLATE(BM_antipode, 4, 9, 1);
BENCHMARK_TEMPLATE(BM_antipode, 4, 9, 2);
BENCHMARK_TEMPLATE(BM_antipode, 4, 9, 3);
//BENCHMARK_TEMPLATE(BM_antipode, 4, 9, 4);

BENCHMARK_TEMPLATE(BM_antipode, 4, 10, 1);
BENCHMARK_TEMPLATE(BM_antipode, 4, 10, 2);
BENCHMARK_TEMPLATE(BM_antipode, 4, 10, 3);
//BENCHMARK_TEMPLATE(BM_antipode, 4, 10, 4);

BENCHMARK_TEMPLATE(BM_antipode, 4, 11, 1);
BENCHMARK_TEMPLATE(BM_antipode, 4, 11, 2);
BENCHMARK_TEMPLATE(BM_antipode, 4, 11, 3);
//BENCHMARK_TEMPLATE(BM_antipode, 4, 11, 4);
//BENCHMARK_TEMPLATE(BM_antipode, 4, 11, 5);

BENCHMARK_TEMPLATE(BM_antipode, 4, 12, 1);
BENCHMARK_TEMPLATE(BM_antipode, 4, 12, 2);
BENCHMARK_TEMPLATE(BM_antipode, 4, 12, 3);
//BENCHMARK_TEMPLATE(BM_antipode, 4, 12, 4);
//BENCHMARK_TEMPLATE(BM_antipode, 4, 12, 5);

BENCHMARK_TEMPLATE(BM_antipode, 4, 13, 1);
BENCHMARK_TEMPLATE(BM_antipode, 4, 13, 2);
BENCHMARK_TEMPLATE(BM_antipode, 4, 13, 3);
//BENCHMARK_TEMPLATE(BM_antipode, 4, 13, 4);
//BENCHMARK_TEMPLATE(BM_antipode, 4, 13, 5);
//BENCHMARK_TEMPLATE(BM_antipode, 4, 13, 6);
