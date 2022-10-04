#include <benchmark/benchmark.h>
#include <iostream>

#include <libalgebra/libalgebra.h>
#include <libalgebra/giles_multiplication.h>
#include <libalgebra/multiplication_helpers.h>



#include <libalgebra/tensor.h>
#include "../../tests/common/random_vector_generator.h"
#include "../../tests/common/rng.h"
#include "../../tests/common/SHOW.h"

using alg::LET;

typedef alg::coefficients::coefficient_field<float> float_field;
typedef alg::coefficients::coefficient_field<alg::coefficients::rational> rational_field;

using namespace alg;

template <DEG CalcDepth, typename Coeffs, DEG Width, DEG Depth, DEG TileLetters, typename Fn>
void fma_traditional(dtl::GilesMultiplier<Coeffs,Width, Depth, TileLetters>& helper, Fn op) noexcept
{
    const auto lhs_deg = static_cast<int>(helper.lhs_degree());
    const auto rhs_deg = static_cast<int>(helper.rhs_degree());

    for (int out_deg = int(CalcDepth); out_deg >= 0; --out_deg) {
        int lhs_deg_max = std::min(out_deg, lhs_deg);
        int lhs_deg_min = std::max(0, out_deg - rhs_deg);

        auto* out_ptr = helper.fwd_write(out_deg);

        for (int lh_deg = lhs_deg_max; lh_deg >= lhs_deg_min; --lh_deg) {
            int rh_deg = out_deg - lh_deg;

            const auto* lhs_ptr = helper.left_fwd_read(static_cast<DEG>(lh_deg), 0);
            const auto* rhs_ptr = helper.right_fwd_read(static_cast<DEG>(rh_deg), 0);

            const auto lhs_size = static_cast<std::ptrdiff_t>(helper.powers[lh_deg]);
            const auto rhs_size = static_cast<std::ptrdiff_t>(helper.powers[rh_deg]);

            auto* ptr = out_ptr;
            for (std::ptrdiff_t i=0; i<lhs_size; ++i) {
                for (std::ptrdiff_t j=0; j<rhs_size; ++j) {
                    *(ptr++) += op(lhs_ptr[i]*rhs_ptr[j]);
                }
            }
        }
    }
}

template <typename Coeffs, DEG Width, DEG Depth>
using dense_ft_vec_t = alg::free_tensor<Coeffs, Width, Depth, alg::vectors::dense_vector>;

template <DEG TileLetters, typename Coeffs, DEG Width, DEG Depth, typename Fn>
void tiled_fma(dense_ft_vec_t<Coeffs, Width, Depth>& out,
               const dense_ft_vec_t<Coeffs, Width, Depth>& lhs_in,
               const dense_ft_vec_t<Coeffs, Width, Depth>& rhs_in,
               Fn op)
{
    static_assert(2*TileLetters <= Depth, "letters in tile cannot exceed half the depth");

    dtl::GilesMultiplier<Coeffs, Width, Depth, TileLetters> helper(out, lhs_in, rhs_in);

    const auto max_depth = static_cast<int>(Depth);
    auto* tile = helper.out_tile_ptr();
    const auto* left_rtile = helper.left_read_tile_ptr();
    const auto* right_rtile = helper.right_read_tile_ptr();

    for (DEG out_deg = max_depth; out_deg > 2 * TileLetters; --out_deg) {
        const DIMN stride = helper.powers[out_deg - TileLetters];

        for (DIMN k = 0; k < helper.powers[out_deg - 2 * TileLetters]; ++k) {
            auto k_reverse = helper.reverse_key(out_deg, k);

            // Handle 0*out_depth and out_depth*0
            const auto& lhs_unit = helper.left_unit();
            const auto* rhs_ptr = helper.right_fwd_read(out_deg, k);
            for (DIMN i = 0; i < helper.tile_width; ++i) {
                for (DIMN j = 0; j < helper.tile_width; ++j) {
                    tile[i * helper.tile_width + j] = op(lhs_unit * rhs_ptr[i * stride + j]);
                }
            }

            const auto& rhs_unit = helper.right_unit();
            const auto* lhs_ptr = helper.left_fwd_read(out_deg, k);
            for (DIMN i = 0; i < helper.tile_width; ++i) {
                for (DIMN j = 0; j < helper.tile_width; ++j) {
                    tile[i * helper.tile_width + j] += op(lhs_ptr[i * stride + j] * rhs_unit);
                }
            }

            // Handle left hand too small cases
            // We iterate through lh_deg < TileLetters, which are the number of
            // letters that we read from the left.
            // (llk,lrk,k,t)
            // llk = left part of left key of degree lh_deg
            // lrk = right part of left key of degree TileLetters - lh_deg
            // k = indexing key of degree out_deg - 2*TileLetters
            // t=  tile key of degree TileLetters
            // lh_deg + (TileLetters-lh_deg) + (out_deg - 2*TileLetters) + TileLetters
            //      = out_deg (good).
            for (DEG lh_deg = 1; lh_deg < TileLetters; ++lh_deg) {
                auto rh_deg = out_deg - (2 * TileLetters + lh_deg);
                for (DIMN i = 0; i < helper.tile_width; ++i) {
                    const auto split = helper.split_key(lh_deg, i);
                    const auto& left_val = *helper.left_fwd_read(lh_deg, split.first);
                    helper.read_right_tile(out_deg - rh_deg, helper.combine_keys(out_deg - 2 * TileLetters, split.second, k));

                    for (DIMN j = 0; j < helper.tile_width; ++j) {
                        tile[i * helper.tile_width + j] += op(left_val * right_rtile[j]);
                    }
                }
            }

            // Handle middle cases
            for (DEG lhs_deg = 0; lhs_deg <= out_deg - 2 * TileLetters; ++lhs_deg) {
                const auto rhs_deg = out_deg - 2 * TileLetters - lhs_deg;
                auto split = helper.split_key(rhs_deg, k);
                //                    assert(split.first*integer_maths::power(WIDTH_HALL_BASIS, rhs_deg) + split.second == k);
                helper.read_left_tile(lhs_deg + TileLetters, helper.reverse_key(lhs_deg, split.first));
                helper.read_right_tile(rhs_deg + TileLetters, split.second);

                for (DIMN i = 0; i < helper.tile_width; ++i) {
                    for (DIMN j = 0; j < helper.tile_width; ++j) {
                        tile[i * helper.tile_width + j] += op(left_rtile[i] * right_rtile[j]);
                    }
                }
            }

            // Handle right hand too small cases
            // (t,k,rlk,rrk)
            // t = tile key of degree TileLetters
            // k = indexing key of degree out_deg - 2*TileLetters
            // rlk = left part of right key of degree TileLetters - rh_deg
            // rrk = right part of right key of degree rh_deg
            for (DEG rh_deg = 1; rh_deg < TileLetters; ++rh_deg) {
                auto lh_deg = out_deg - (2 * TileLetters + rh_deg);
                for (DIMN j = 0; j < helper.tile_width; ++j) {
                    const auto split = helper.split_key(rh_deg, j);
                    const auto& right_val = *helper.right_fwd_read(rh_deg, split.second);
                    helper.read_left_tile(out_deg - rh_deg, helper.combine_keys(TileLetters - rh_deg, k_reverse, helper.reverse_key(TileLetters - rh_deg, split.first)));

                    for (DIMN i = 0; i < helper.tile_width; ++i) {
                        tile[helper.reverse(i) * helper.tile_width + j] += op(left_rtile[i] * right_val);
                    }
                }
            }

            // Finally, write out the answer to the output buffers
            helper.write_tile(out_deg, k, k_reverse);
        }
    }

    fma_traditional<2*TileLetters>(helper, op);
}


template<unsigned Width, unsigned Depth>
class RandomRationalDenseFixture : public ::benchmark::Fixture {
public:

    static constexpr unsigned width = Width;
    static constexpr unsigned depth = Depth;

    typedef alg::free_tensor<rational_field, Width, Depth, alg::vectors::dense_vector> TENSOR;
    typedef alg::lie<rational_field, Width, Depth, alg::vectors::dense_vector> LIE;

    using rat_dist = la_testing::uniform_rational_distribution<rational_field::S>;
    using rvg_t = la_testing::random_vector_generator<TENSOR, rat_dist>;
    using rvg_l = la_testing::random_vector_generator<LIE, rat_dist>;

    const TENSOR tunit;
    const TENSOR tzero;

    std::mt19937 rngt;
    rvg_t rvgt;

    std::mt19937 rngl;
    rvg_l rvgl;

    typedef typename TENSOR::KEY KEY;

    RandomRationalDenseFixture() : tunit(KEY()), tzero(),
                                   rngt(std::random_device()()), rvgt(-1, 1),
                                   rngl(std::random_device()()), rvgl(-1, 1)

    {}

};


BENCHMARK_TEMPLATE_DEFINE_F(RandomRationalDenseFixture, TraditionalMultiplication, benchmark::State.range(0), benchmark::State.range(1))(benchmark::State& state) {

    auto lhs = rvgt(rngt);
    auto rhs = rvgt(rngt);

    const DEG TileLetters = 1;        // TODO: implement properly

    for (auto _ : state) {
        auto result = lhs*rhs;
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}

BENCHMARK_TEMPLATE_DEFINE_F(RandomRationalDenseFixture, TiledMultiplication, 4, 4)(benchmark::State& state) {

    auto lhs = rvgt(rngt);
    auto rhs = rvgt(rngt);

    //    SHOW(lhs);
    //    SHOW(rhs);

    TENSOR result;

    const DEG TileLetters = 1;        // TODO: implement properly

    dtl::GilesMultiplier<rational_field , 4, 4, TileLetters> helper(result, lhs, rhs);

    for (auto _ : state) {
        tiled_fma<1>(result, lhs, rhs, mult::scalar_passthrough());
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}

static void CustomArguments(benchmark::internal::Benchmark* b) {
    // define I, J
//    for (auto i : I)
    for (int i = 5; i <= 10; i++)
    {
//        for (auto j : J)
    for (int j = 5; j <= 10; j++)
            b->Args({i, j});
    }
}

BENCHMARK_REGISTER_F(RandomRationalDenseFixture, TraditionalMultiplication)->Apply(CustomArguments);

BENCHMARK_REGISTER_F(RandomRationalDenseFixture, TiledMultiplication);//->Apply(CustomArguments);

BENCHMARK_MAIN();