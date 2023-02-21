#include <UnitTest++.h>

#include <iostream>
#include <vector>

#include <libalgebra/libalgebra.h>
#include <libalgebra/tensor.h>
#include "../../common/random_vector_generator.h"
#include "../../common/rng.h"

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

//static SHUFFLE_TENSOR shift_down(const SHUFFLE_TENSOR& sh, typename SHUFFLE_TENSOR::KEY word)
//{
//    SHUFFLE_TENSOR result(sh), working;
//    while (word.size()) {
//        auto letter = word.lparent();
//        word = word.rparent();
//        for (auto pr : result) {
//            if (pr.key().size() > 0 && pr.key().lparent() == letter)
//                working[pr.key().rparent()] = result[pr.key()];
//        }
//        result.swap(working);
//        working.clear();
//    }
//    return result;
//}
//
//// the evaluation of the adjoint operation is worked out below
//// <sh,ab>=\sum_{uv=sh}<ua><vb>
////        = <\sum_{uv=sh}<ua>v,b>
//// Let T_w(sh) be all the projection of sh onto the part beginning with w with w removed
//// \sum_{i} <ki shi,ab> =
////        = \sum_{i} ki<\sum_{uv=shi}<ua>v,b>
////        = \sum_{u} < <ua>T_u(sh), b>
////  The action of the adjoint of tensor multiplication by a multiplication is \sum_u <ua> T_u(sh)
//
//static SHUFFLE_TENSOR adjoint_to_multiply(const TENSOR& t, SHUFFLE_TENSOR sh)
//{
//    // this implementation is understandable and reliable but repetitive and can be radically accelerated
//    SHUFFLE_TENSOR result;
//    for (auto& pr : t) {
//        result += shift_down(sh, pr.key()) * pr.value();
//    }
//    return result;
//}

SUITE(Adjoint)
{
    typedef DenseFixture<float_field, 4, 4> dense_fixture;
    TEST_FIXTURE(dense_fixture, Test)
    {
        CHECK_EQUAL(0,1);
    }

}// SUITE adjoint
