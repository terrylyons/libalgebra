#include <UnitTest++/UnitTest++.h>

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

SUITE(Multiplication)
{
    typedef DenseFixture<float_field, 4, 4> dense_fixture;

    TEST_FIXTURE(dense_fixture, cmake_check)
    {
        CHECK_EQUAL(0, 1);
    }

}// SUITE Multiplication