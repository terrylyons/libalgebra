#include <UnitTest++/UnitTest++.h>

#include <iostream>
#include <vector>
#include "SHOW.h"
#include <functional>
#include <sstream>
#include <string>

#include <libalgebra/libalgebra.h>
#include <libalgebra/tensor.h>
#include "../../common/random_vector_generator.h"
#include "../../common/rng.h"

#include <libalgebra/giles_multiplication.h>

using alg::LET;

typedef alg::coefficients::coefficient_field<float> float_field;
typedef alg::coefficients::coefficient_field<alg::coefficients::rational> rational_field;

template<typename Coeff, unsigned Width, unsigned Depth>
struct SparseFixture {

    typedef typename Coeff::S S;
    typedef typename Coeff::Q Q;

    typedef alg::free_tensor_basis<Width, Depth> TBASIS;
    typedef alg::vectors::sparse_vector<TBASIS, Coeff> VECT;
    typedef alg::free_tensor<Coeff, Width, Depth, alg::vectors::sparse_vector> TENSOR;

    typedef typename TBASIS::KEY KEY;

    const TENSOR tunit;
    const TENSOR tzero;

    SparseFixture() : tunit(KEY()), tzero()
    {}

    KEY make_key(const LET* arg, const std::size_t N)
    {
        KEY k;
        for (std::size_t i = 0; i < N; ++i) {
            k.push_back(arg[i]);
        }
        return k;
    }
};

template<typename Coeff, unsigned Width, unsigned Depth>
struct DenseFixture {

    typedef typename Coeff::S S;
    typedef typename Coeff::Q Q;

    static constexpr unsigned width = Width;
    static constexpr unsigned depth = Depth;
    typedef Coeff coeffs;

    typedef alg::free_tensor_basis<Width, Depth> TBASIS;
    typedef alg::vectors::dense_vector<TBASIS, Coeff> VECT;
    typedef alg::free_tensor<Coeff, Width, Depth, alg::vectors::dense_vector> TENSOR;

    typedef typename TBASIS::KEY KEY;

    const TENSOR tunit;
    const TENSOR tzero;

    DenseFixture() : tunit(KEY()), tzero()
    {}

    KEY make_key(const LET* arg, const std::size_t N)
    {
        KEY k;
        for (std::size_t i = 0; i < N; ++i) {
            k.push_back(arg[i]);
        }
        return k;
    }
};

template<unsigned Width, unsigned Depth>
struct RandomRationalDenseFixture {

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

using namespace alg;

template<DEG WIDTH, DEG DEPTH>
struct Environment {
    using scalar_field = coefficients::rational_field;
    using S = typename scalar_field::S;

    // note name changes - error: declaration of 'constexpr const DEG DEPTH' shadows template parameter

    static constexpr DEG WIDTH_HALL_BASIS = WIDTH; // name change
    static constexpr DEG DEPTH_HALL_BASIS = DEPTH; // name change
    static constexpr DEG poly_width = hall_basis<WIDTH_HALL_BASIS, DEPTH_HALL_BASIS>::start_of_degree(DEPTH + 1);

    using poly_t = alg::poly<scalar_field>;
    using poly_coeffs = coefficients::coefficient_ring<poly_t, typename scalar_field::Q>;
    static_assert(std::is_same<poly_coeffs::S, poly_t>::value, "the trait class of a coefficient ring must report the type of the coefficients");

    template<DEG DEPTH_LIE> // name change
    using LIE_ = alg::lie<poly_coeffs, WIDTH, DEPTH_LIE, vectors::dense_vector>; // name change
    using LIE = LIE_<DEPTH>;

    template<DEG DEPTH_LIE_BASIS> // name change
    using lie_basis_t_ = lie_basis<WIDTH, DEPTH_LIE_BASIS>; // name change

    using lie_basis_t = lie_basis_t_<WIDTH>;
    lie_basis_t lbasis;

    template<DEG DEPTH_TENSOR> // name change
    using TENSOR_ = alg::free_tensor<poly_coeffs, WIDTH, DEPTH_TENSOR, vectors::dense_vector>; // name change
    using TENSOR = TENSOR_<DEPTH>;

    template<DEG DEPTH_TENSOR_BASIS> // name change
    using tensor_basis_t_ = alg::tensor_basis<WIDTH, DEPTH_TENSOR_BASIS>; // name change

    using tensor_basis_t = tensor_basis_t_< DEPTH>;
    tensor_basis_t tbasis;

    //using SHUFFLE_TENSOR = alg::shuffle_tensor<scalar_field, WIDTH, DEPTH>;
    template<DEG DEPTH_SHUFFLE_TENSOR> // name change
    using SHUFFLE_TENSOR_ = alg::shuffle_tensor<poly_coeffs, WIDTH, DEPTH_SHUFFLE_TENSOR>; // name change
    using SHUFFLE_TENSOR = SHUFFLE_TENSOR_<DEPTH>;

    template<DEG DEPTH_SHUFFLE_TENSOR_BASIS> // name change
    using shuffle_tensor_basis_t_ = alg::tensor_basis<WIDTH, DEPTH_SHUFFLE_TENSOR_BASIS>; // name change
    using shuffle_tensor_basis_t = shuffle_tensor_basis_t_<DEPTH>;
    shuffle_tensor_basis_t sbasis;

    using SHUFFLE_TENSOR_OVER_POLYS = alg::shuffle_tensor<poly_coeffs, WIDTH, DEPTH>;
    using inner_product = alg::operators::shuffle_tensor_functional<TENSOR, SHUFFLE_TENSOR_OVER_POLYS>;

    using MAPS = maps<poly_coeffs, WIDTH, DEPTH, TENSOR, LIE>;
    using CBH = cbh<poly_coeffs, WIDTH, DEPTH, TENSOR, LIE>;

    MAPS maps_;
    CBH cbh_;

    // make a LIE element whose hall coefficients are monomials
    LIE generic_lie(const int offset = 0) const
    {
        LIE result;
        for (auto lie_key : lbasis.iterate_keys()) {
            result.add_scal_prod(lie_key, poly_t(lie_key + offset, S(1)));
        }
        return result;
    }

    // The bilinear function K takes a shuffle and contracts it with a tensor to produce a scalar.
    // In this case the scalar is itself a polynomial
    static typename poly_coeffs::S K(const SHUFFLE_TENSOR& functional, const TENSOR& sig_data)
    {
        SHUFFLE_TENSOR_OVER_POLYS functional_p;
        for (auto& key_value : functional)
            functional_p.add_scal_prod(key_value.key(), typename poly_coeffs::S(key_value.value()));
        inner_product f(sig_data);//todo: left and right switched here?
        return f(functional_p);
    }

    // creates a generic vector with monomial coefficients
    template<class VECTOR, const int offset0 = 0>
    static VECTOR generic_vector(const int offset = offset0)
    {
        const typename VECTOR::BASIS& basis = VECTOR::basis;
        // LIE basis starts at 1 which is confusing

        int count;
        if constexpr (std::is_integral<decltype(basis.begin())>::value) {
            // invalid code if the premise is false - constexpr essential to avoid compilation
            count = basis.begin();
        }
        else {
            count = 0;
        }

        std::map<int, std::pair<typename VECTOR::KEY, std::string>> legend;

        VECTOR result;
        //for range type for loop approach use ": basis.iterate_keys()"
        for (auto key = basis.begin(), end = basis.end(); key != end; key = basis.nextkey(key)) {
            result[key] = poly_t(count + offset, 1);

            // record the mapping from keys to monomials
            auto basis_key_pair = std::pair<typename VECTOR::BASIS*, typename VECTOR::KEY>(&VECTOR::basis, key);
            std::stringstream buffer;
            buffer << basis_key_pair;
//            legend[count + offset] = std::pair(key, buffer.str());  // TODO: error: missing template arguments before '(' token
            legend[count + offset] = std::pair<typename VECTOR::KEY, std::string>(key, buffer.str());
            std::cout << " monomial index:" << count + offset << " basis index:" << key << " basis value:" << buffer.str() << "\n";

            ++count;
        }
        return result;
    }


};

template<class FULL_TENSOR, class SHORT_TENSOR>
FULL_TENSOR& add_equals_short(FULL_TENSOR& out, const SHORT_TENSOR& in)
{
    const static typename FULL_TENSOR::BASIS fbasis;
    const static typename SHORT_TENSOR::BASIS sbasis;
    auto key = sbasis.begin(), end = sbasis.end();
    auto fkey = fbasis.begin(), fend = fbasis.end();
    for (; key != end && fkey != fend; key = sbasis.nextkey(key), fkey = fbasis.nextkey(fkey))
        out[fkey] += in[key];
    return out;
}

template<class EnvironmentOUT, class EnvironmentIN>
typename EnvironmentOUT::TENSOR apply1(const std::map<typename EnvironmentOUT::TENSOR::BASIS::KEY, typename EnvironmentIN::SHUFFLE_TENSOR>& result, const typename EnvironmentIN::TENSOR& in)
{
    typename EnvironmentOUT::TENSOR out;
    for (const auto& x : result) {
        auto& key = x.first;
        auto& tvalue = x.second;
        // both environments must have compatible scalar types
        out[key] = typename EnvironmentIN::K(tvalue, in);
    };

    return out;
}

SUITE(Multiplication)
{
    typedef DenseFixture<float_field, 4, 4> dense_fixture;

    TEST_FIXTURE(dense_fixture, giles_multiplication_check)
    {
        LET k1[] = {1,2};

        TENSOR lhs(make_key(k1, 2));
        TENSOR rhs(make_key(k1, 2));

        vectors::dtl::vector_base_access::convert(lhs).resize_to_dimension(341); // TODO: implement properly
        vectors::dtl::vector_base_access::convert(rhs).resize_to_dimension(341); // TODO: implement properly

        TENSOR expected = lhs*rhs;

        TENSOR result;

        const DEG TileLetters = 1;        // TODO: implement properly
        DEG tensor_width = width;        // TODO: implement properly
        DEG max_depth = depth;         // TODO: implement properly

        dtl::GilesMultiplier<float_field , width, depth, TileLetters> helper(result, lhs, rhs);

        // objects that come from GilesMultiplier: helper.powers, helper.reverse_key, helper.tile_width,
        // helper.left_forward_read_ptr, helper.reverse, helper.split, tile, left_rtile, right_rtile

        auto* tile = helper.out_tile_ptr();
        const auto* left_rtile = helper.left_read_tile_ptr();
        const auto* right_rtile = helper.right_read_tile_ptr();

        // ############## from multiplication.cpp ################## //


        for (DEG out_deg = max_depth; out_deg > 2 * TileLetters; --out_deg) {
            const DIMN stride = integer_maths::power(tensor_width, out_deg - TileLetters);

            for (DIMN k = 0; k < helper.powers[out_deg - 2 * TileLetters]; ++k) {
                auto k_reverse = helper.reverse_key(out_deg, k);

                // Handle 0*out_depth and out_depth*0
                const auto& lhs_unit = helper.left_unit();
                const auto* rhs_ptr = helper.right_fwd_read(out_deg, k);
                for (DIMN i = 0; i < helper.tile_width; ++i) {
                    for (DIMN j = 0; j < helper.tile_width; ++j) {
                        tile[i * helper.tile_width + j] = lhs_unit * rhs_ptr[i * stride + j];
                    }
                }

                const auto& rhs_unit = helper.right_unit();
                const auto* lhs_ptr = helper.left_fwd_read(out_deg, k);
                for (DIMN i = 0; i < helper.tile_width; ++i) {
                    for (DIMN j = 0; j < helper.tile_width; ++j) {
                        tile[i * helper.tile_width + j] += lhs_ptr[i * stride + j] * rhs_unit;
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
                    for (DIMN i = 0; i < helper.powers[lh_deg]; ++i) {
                        const auto split = helper.split_key(lh_deg, i);
                        const auto& left_val = *helper.left_fwd_read(lh_deg, split.first);
                        helper.read_right_tile(out_deg - rh_deg, helper.combine_keys(out_deg - 2 * TileLetters, split.second, k));

                        for (DIMN j = 0; j < helper.tile_width; ++j) {
                            tile[i * helper.tile_width + j] += left_val * right_rtile[j];
                        }
                    }
                }

                // Handle middle cases
                for (DEG lhs_deg = 0; lhs_deg <= out_deg - 2 * TileLetters; ++lhs_deg) {
                    const auto rhs_deg = out_deg - 2 * TileLetters - lhs_deg;
                    auto split = helper.split_key(lhs_deg, k);
                    helper.read_left_tile(lhs_deg + TileLetters, helper.reverse_key(lhs_deg, split.first));
                    helper.read_right_tile(rhs_deg + TileLetters, split.second);

                    for (DIMN i = 0; i < helper.tile_width; ++i) {
                        for (DIMN j = 0; j < helper.tile_width; ++j) {
                            tile[helper.reverse(i) * helper.tile_width + j] += left_rtile[i] * right_rtile[j];
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
                    for (DIMN j = 0; j < helper.powers[rh_deg]; ++j) {
                        const auto split = helper.split_key(rh_deg, j);
                        const auto& right_val = *helper.right_fwd_read(rh_deg, split.second);
                        helper.read_left_tile(out_deg - rh_deg, helper.combine_keys(TileLetters - rh_deg, k_reverse, helper.reverse_key(TileLetters - rh_deg, split.first)));

                        for (DIMN i = 0; i < helper.tile_width; ++i) {
                            tile[helper.reverse(i) * helper.tile_width + j] += left_rtile[i] * right_val;
                        }
                    }
                }

                // Finally, write out the answer to the output buffers
                helper.write_tile(out_deg, k, k_reverse);
            }
        }
        std::cout << "lhs=" << lhs << std::endl;
        std::cout << "rhs=" << rhs << std::endl;
        std::cout << "result=" << result << std::endl;

        // ############## end multiplication.cpp ################## //

        CHECK_EQUAL(expected, result);

    }

    using IN = Environment<4, 4>;
    using SHORT_LIE = IN::LIE_<4>;

    TEST_FIXTURE(IN, giles_multiplication_check_with_polynomials)
    {

        IN in;
        const auto& sbasis = in.sbasis;

        IN::LIE logsig_lhs;
        add_equals_short(logsig_lhs, in.generic_vector<IN::LIE>(1000));
//        SHOW(logsig_lhs);
        IN::TENSOR tensor_logsig_lhs = in.maps_.l2t(logsig_lhs);
        IN::TENSOR lhs = exp(tensor_logsig_lhs);
//        SHOW(lhs);

        IN::LIE logsig_rhs;
        add_equals_short(logsig_rhs, in.generic_vector<IN::LIE>(1000));
//        SHOW(logsig_rhs);
        IN::TENSOR tensor_logsig_rhs = in.maps_.l2t(logsig_rhs);
        IN::TENSOR rhs = exp(tensor_logsig_rhs);
//        SHOW(rhs);

        TENSOR expected = lhs*rhs;

        TENSOR result;

        vectors::dtl::vector_base_access::convert(lhs).resize_to_dimension(341); // TODO: implement properly
        vectors::dtl::vector_base_access::convert(rhs).resize_to_dimension(341); // TODO: implement properly

        const DEG TileLetters = 1;        // TODO: implement properly
        DEG tensor_width = 4;        // TODO: implement properly
        DEG max_depth = 4;         // TODO: implement properly

        dtl::GilesMultiplier<poly_coeffs , 4, 4, TileLetters> helper(result, lhs, rhs);

        // objects that come from GilesMultiplier: helper.powers, helper.reverse_key, helper.tile_width,
        // helper.left_forward_read_ptr, helper.reverse, helper.split, tile, left_rtile, right_rtile

        auto* tile = helper.out_tile_ptr();
        const auto* left_rtile = helper.left_read_tile_ptr();
        const auto* right_rtile = helper.right_read_tile_ptr();

        // ############## from multiplication.cpp ################## //


        for (DEG out_deg = max_depth; out_deg > 2 * TileLetters; --out_deg) {
            const DIMN stride = integer_maths::power(tensor_width, out_deg - TileLetters);

            for (DIMN k = 0; k < helper.powers[out_deg - 2 * TileLetters]; ++k) {
                auto k_reverse = helper.reverse_key(out_deg, k);

                // Handle 0*out_depth and out_depth*0
                const auto& lhs_unit = helper.left_unit();
                const auto* rhs_ptr = helper.right_fwd_read(out_deg, k);
                for (DIMN i = 0; i < helper.tile_width; ++i) {
                    for (DIMN j = 0; j < helper.tile_width; ++j) {
                        tile[i * helper.tile_width + j] = lhs_unit * rhs_ptr[i * stride + j];
                    }
                }

                const auto& rhs_unit = helper.right_unit();
                const auto* lhs_ptr = helper.left_fwd_read(out_deg, k);
                for (DIMN i = 0; i < helper.tile_width; ++i) {
                    for (DIMN j = 0; j < helper.tile_width; ++j) {
                        tile[i * helper.tile_width + j] += lhs_ptr[i * stride + j] * rhs_unit;
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
                    for (DIMN i = 0; i < helper.powers[lh_deg]; ++i) {
                        const auto split = helper.split_key(lh_deg, i);
                        const auto& left_val = *helper.left_fwd_read(lh_deg, split.first);
                        helper.read_right_tile(out_deg - rh_deg, helper.combine_keys(out_deg - 2 * TileLetters, split.second, k));

                        for (DIMN j = 0; j < helper.tile_width; ++j) {
                            tile[i * helper.tile_width + j] += left_val * right_rtile[j];
                        }
                    }
                }

                // Handle middle cases
                for (DEG lhs_deg = 0; lhs_deg <= out_deg - 2 * TileLetters; ++lhs_deg) {
                    const auto rhs_deg = out_deg - 2 * TileLetters - lhs_deg;
                    auto split = helper.split_key(lhs_deg, k);
                    helper.read_left_tile(lhs_deg + TileLetters, helper.reverse_key(lhs_deg, split.first));
                    helper.read_right_tile(rhs_deg + TileLetters, split.second);

                    for (DIMN i = 0; i < helper.tile_width; ++i) {
                        for (DIMN j = 0; j < helper.tile_width; ++j) {
                            tile[helper.reverse(i) * helper.tile_width + j] += left_rtile[i] * right_rtile[j];
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
                    for (DIMN j = 0; j < helper.powers[rh_deg]; ++j) {
                        const auto split = helper.split_key(rh_deg, j);
                        const auto& right_val = *helper.right_fwd_read(rh_deg, split.second);
                        helper.read_left_tile(out_deg - rh_deg, helper.combine_keys(TileLetters - rh_deg, k_reverse, helper.reverse_key(TileLetters - rh_deg, split.first)));

                        for (DIMN i = 0; i < helper.tile_width; ++i) {
                            tile[helper.reverse(i) * helper.tile_width + j] += left_rtile[i] * right_val;
                        }
                    }
                }

                // Finally, write out the answer to the output buffers
                helper.write_tile(out_deg, k, k_reverse);
            }
        }
        std::cout << "lhs=" << lhs << std::endl;
        std::cout << "rhs=" << rhs << std::endl;
        std::cout << "result=" << result << std::endl;

        // ############## end multiplication.cpp ################## //

        CHECK_EQUAL(expected, result);
    }

}// SUITE Multiplication