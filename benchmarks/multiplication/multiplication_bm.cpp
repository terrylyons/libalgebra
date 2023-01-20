#include <benchmark/benchmark.h>
#include <iostream>

#include "../../tests/common/random_vector_generator.h"
#include "../../tests/common/rng.h"
#include <libalgebra/libalgebra.h>
#include <libalgebra/tensor.h>

using namespace alg;
using alg::LET;
typedef alg::coefficients::coefficient_field<float> float_field;
typedef alg::coefficients::coefficient_field<alg::coefficients::rational> rational_field;
using float_dist = std::uniform_real_distribution<float_field::S>;

template<DEG Width, DEG Depth>
static void BM_traditional_multiplication(benchmark::State& state)
{

    typedef alg::free_tensor<float_field, Width, Depth, alg::vectors::dense_vector> TENSOR;
    typedef alg::lie<float_field, Width, Depth, alg::vectors::dense_vector> LIE;

    typedef typename TENSOR::KEY KEY;

    TENSOR lhs, rhs;
    lhs.base_vector().resize_to_degree(Depth);
    rhs.base_vector().resize_to_degree(Depth);

    KEY kunit;
    lhs[kunit] = 0.983453f;
    rhs[kunit] = 1.0023524f;

    using mul_t = traditional_free_tensor_multiplication<Width, Depth>;
    using traits = dtl::multiplication_traits<mul_t>;
    mul_t mul;

    lhs.base_vector().construct_reverse_data(Depth - 1);// Optional

    TENSOR result;

    for (auto _ : state) {
        benchmark::DoNotOptimize(result.base_vector().as_mut_ptr());
        traits::multiply_and_add(mul, result, lhs, rhs);
        benchmark::ClobberMemory();
    }

    state.SetBytesProcessed(3 * sizeof(float) * state.iterations() * lhs.base_vector().dimension());
}

template<DEG Width, DEG Depth, DEG TileLetters, DEG CacheLetters=0>
static void BM_tiled_multiplication(benchmark::State& state)
{

    //    TODO: replace random vector generator for dense vectors
    //    auto* ptr = lhs.base_vector().as_mut_ptr();
    //    std::generate(ptr, ptr + lhs.base_vector().dimension(), [this, distribution]() { return distribution(rng);

    typedef alg::free_tensor<float_field, Width, Depth, alg::vectors::dense_vector> TENSOR;
    typedef alg::lie<float_field, Width, Depth, alg::vectors::dense_vector> LIE;

    typedef typename TENSOR::KEY KEY;

    TENSOR lhs, rhs;
    lhs.base_vector().resize_to_degree(Depth);
    rhs.base_vector().resize_to_degree(Depth);

    KEY kunit;
    lhs[kunit] = 0.983453f;
    rhs[kunit] = 1.0023524f;

    using mul_t = tiled_free_tensor_multiplication<Width, Depth, TileLetters, CacheLetters>;
    using traits = dtl::multiplication_traits<mul_t>;
    mul_t mul;

    lhs.base_vector().construct_reverse_data(Depth - 1);// Optional

    TENSOR result;

    for (auto _ : state) {
        benchmark::DoNotOptimize(result.base_vector().as_mut_ptr());
        traits::multiply_and_add(mul, result, lhs, rhs);
        benchmark::ClobberMemory();
    }

    state.SetBytesProcessed(3 * sizeof(float) * state.iterations() * lhs.base_vector().dimension());
}

//BENCHMARK_TEMPLATE(BM_traditional_multiplication, 4, 4);
//BENCHMARK_TEMPLATE(BM_traditional_multiplication, 4, 5);
//BENCHMARK_TEMPLATE(BM_traditional_multiplication, 4, 6);
//BENCHMARK_TEMPLATE(BM_traditional_multiplication, 4, 7);
//BENCHMARK_TEMPLATE(BM_traditional_multiplication, 4, 8);
//BENCHMARK_TEMPLATE(BM_traditional_multiplication, 4, 9);
//BENCHMARK_TEMPLATE(BM_traditional_multiplication, 4, 10);
//BENCHMARK_TEMPLATE(BM_traditional_multiplication, 4, 11);
//BENCHMARK_TEMPLATE(BM_traditional_multiplication, 4, 12);
BENCHMARK_TEMPLATE(BM_traditional_multiplication, 4, 13);

//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 4, 1, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 5, 1, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 6, 1, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 7, 1, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 8, 1, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 9, 1, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 10, 1, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 11, 1, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 12, 1, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 1, 0);

//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 4, 2, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 5, 2, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 6, 2, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 7, 2, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 8, 2, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 9, 2, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 10, 2, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 11, 2, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 12, 2, 0);
BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 2, 0);

//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 4, 3, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 5, 3, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 6, 3, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 7, 3, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 8, 3, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 9, 3, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 10, 3, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 11, 3, 0);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 12, 3, 0);
BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 3, 0);


//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 4, 1, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 5, 1, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 6, 1, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 7, 1, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 8, 1, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 9, 1, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 10, 1, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 11, 1, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 12, 1, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 1, 1);

//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 4, 2, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 5, 2, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 6, 2, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 7, 2, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 8, 2, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 9, 2, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 10, 2, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 11, 2, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 12, 2, 1);
BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 2, 1);

//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 4, 3, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 5, 3, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 6, 3, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 7, 3, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 8, 3, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 9, 3, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 10, 3, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 11, 3, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 12, 3, 1);
BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 3, 1);

//
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 4, 1, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 5, 1, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 6, 1, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 7, 1, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 8, 1, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 9, 1, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 10, 1, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 11, 1, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 12, 1, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 1, 2);

//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 4, 2, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 5, 2, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 6, 2, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 7, 2, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 8, 2, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 9, 2, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 10, 2, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 11, 2, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 12, 2, 2);
BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 2, 2);

//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 4, 3, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 5, 3, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 6, 3, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 7, 3, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 8, 3, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 9, 3, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 10, 3, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 11, 3, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 12, 3, 2);
BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 3, 2);


BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 3, 3);
BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 4, 0);
BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 4, 1);
BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 4, 2);


//
//BENCHMARK_TEMPLATE(BM_traditional_multiplication, 4, 4);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 4, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 4, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 4, 3);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 4, 4);
//
//
//BENCHMARK_TEMPLATE(BM_traditional_multiplication, 4, 5);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 5, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 5, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 5, 3);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 5, 4);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 5, 5);
//
//BENCHMARK_TEMPLATE(BM_traditional_multiplication, 4, 6);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 6, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 6, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 6, 3);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 6, 4);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 6, 5);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 6, 6);
//
//BENCHMARK_TEMPLATE(BM_traditional_multiplication, 4, 7);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 7, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 7, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 7, 3);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 7, 4);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 7, 5);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 7, 6);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 7, 7);
//
//BENCHMARK_TEMPLATE(BM_traditional_multiplication, 4, 8);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 8, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 8, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 8, 3);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 8, 4);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 8, 5);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 8, 6);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 8, 7);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 8, 8);
//
//BENCHMARK_TEMPLATE(BM_traditional_multiplication, 4, 9);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 9, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 9, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 9, 3);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 9, 4);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 9, 5);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 9, 6);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 9, 7);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 9, 8);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 9, 9);
//
//BENCHMARK_TEMPLATE(BM_traditional_multiplication, 4, 10);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 10, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 10, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 10, 3);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 10, 4);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 10, 5);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 10, 6);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 10, 7);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 10, 8);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 10, 9);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 10, 10);
//
//BENCHMARK_TEMPLATE(BM_traditional_multiplication, 4, 11);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 11, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 11, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 11, 3);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 11, 4);
////BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 11, 5); // SIGSEGV
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 11, 6);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 11, 7);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 11, 8);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 11, 9);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 11, 10);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 11, 11);
//
//BENCHMARK_TEMPLATE(BM_traditional_multiplication, 4, 12);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 12, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 12, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 12, 3);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 12, 4);
////BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 12, 5); // SIGSEGV
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 12, 6);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 12, 7);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 12, 8);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 12, 9);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 12, 10);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 12, 11);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 12, 12);
//
//BENCHMARK_TEMPLATE(BM_traditional_multiplication, 4, 13);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 1);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 2);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 3);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 4);
////BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 5); // SIGSEGV
////BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 6); // SIGSEGV
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 7);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 8);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 9);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 10);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 11);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 12);
//BENCHMARK_TEMPLATE(BM_tiled_multiplication, 4, 13, 13);
