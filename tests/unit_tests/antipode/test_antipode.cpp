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

SUITE(Antipode)
{

    // ##################### DENSE TESTS ##################### //

    typedef DenseFixture<float_field, 4, 4> dense_fixture;

    // test: antipode(zero) == zero

    TEST_FIXTURE(dense_fixture, DenseAntipodeZeroTest)
    {
        // check: { } --> { }

        TENSOR result = antipode(tzero);

        CHECK_EQUAL(tzero, result);
    }

    // test: antipode(identity) == identity

    TEST_FIXTURE(dense_fixture, DenseAntipodeIdentiyTest)
    {
        // check: { 1{} } --> { 1{} }

        TENSOR result = antipode(tunit);

        CHECK_EQUAL(tunit, result);
    }

    // test: antipode(length 1 word) == - length 1 word

    TEST_FIXTURE(dense_fixture, DenseAntipodeOneLetterTest)
    {
        // check: { 1{1} } --> { -1{1} }

        LET k1[] = {1};

        TENSOR input_tensor(make_key(k1, 1));

        TENSOR expected;

        expected.add_scal_prod(make_key(k1, 1), -1.0);

        TENSOR result = antipode(input_tensor);

        CHECK_EQUAL(expected, result);
    }

    // test: key/index look-ups for single even word tensor

    TEST_FIXTURE(dense_fixture, DenseAntipodeOneEvenWordTest)
    {

        // check {1{12}} --> {1{21}}

        LET k12[] = {1, 2};
        LET k21[] = {2, 1};

        TENSOR input_tensor(make_key(k12, 2));
        TENSOR expected(make_key(k21, 2));

        TENSOR result = antipode(input_tensor);

        CHECK_EQUAL(expected, result);
    }

    // test: key/index look-ups for single odd word tensor

    TEST_FIXTURE(dense_fixture, DenseAntipodeOneOddWordTest)
    {

        // check {1{123}} --> {-1{321}}

        LET k123[] = {1, 2, 3};
        LET k321[] = {3, 2, 1};

        TENSOR input_tensor(make_key(k123, 3));

        TENSOR expected;
        expected.add_scal_prod(make_key(k321, 3), -1.0);

        TENSOR result = antipode(input_tensor);

        CHECK_EQUAL(expected, result);
    }

    // test: key/index look-ups for multiple word tensor

    TEST_FIXTURE(dense_fixture, DenseAntipodeMultipleWordTest)
    {
            // check: { 1.0{11} 2.0{12} 3.0{21} 4.0{22} 5.0{123} } --> { 1.0{11} 3.0{12} 2.0{21} 4.0{22} -5.0{321} }

        LET k11[] = {1, 1};
        LET k12[] = {1, 2};
        LET k21[] = {2, 1};
        LET k22[] = {2, 2};

        LET k123[] = {1, 2, 3};
        LET k321[] = {3, 2, 1};

        TENSOR input_tensor;

        input_tensor.add_scal_prod(make_key(k11, 2), 1.0);
        input_tensor.add_scal_prod(make_key(k12, 2), 2.0);
        input_tensor.add_scal_prod(make_key(k21, 2), 3.0);
        input_tensor.add_scal_prod(make_key(k22, 2), 4.0);

        input_tensor.add_scal_prod(make_key(k123, 3), 5.0);

        TENSOR expected;

        expected.add_scal_prod(make_key(k11, 2), 1.0);
        expected.add_scal_prod(make_key(k12, 2), 3.0);
        expected.add_scal_prod(make_key(k21, 2), 2.0);
        expected.add_scal_prod(make_key(k22, 2), 4.0);

        expected.add_scal_prod(make_key(k321, 3), -5.0);

        TENSOR result = antipode(input_tensor);

        CHECK_EQUAL(expected, result);
    }

    typedef RandomRationalDenseFixture<4, 4> random_rational_dense_fixture;

    // test using random tensor: antipode(antipode(random_tensor)) == random_tensor

    TEST_FIXTURE(random_rational_dense_fixture, DenseAntipodeTwiceTest)
    {

        auto random_tensor = rvgt(rngt);

        CHECK_EQUAL(random_tensor, antipode(antipode(random_tensor)));

    }

    // tests using group-like elements: antipode(group_like_tensor) * group_like_tensor == Identity

    TEST_FIXTURE(random_rational_dense_fixture, DenseAntipodeGroupLikeIdentityTest)
    {
        auto random_lie = rvgl(rngl);

        alg::maps<rational_field, width, depth, TENSOR, LIE> maps;

        auto group_like_tensor = exp(maps.l2t(random_lie));

        auto result = antipode(group_like_tensor) * group_like_tensor;

        CHECK_EQUAL(tunit, result);

    }

    // ##################### SPARSE TESTS ##################### //

    typedef SparseFixture<float_field, 4, 4> sparse_fixture;

    // test: key/index look-ups for multiple word tensor

    TEST_FIXTURE(sparse_fixture, SparseAntipodeMultipleWordTest)
    {
        // check: { 1.0{11} 2.0{12} 3.0{21} 4.0{22} 5.0{123} } --> { 1.0{11} 3.0{12} 2.0{21} 4.0{22} -5.0{321} }

        LET k11[] = {1, 1};
        LET k12[] = {1, 2};
        LET k21[] = {2, 1};
        LET k22[] = {2, 2};

        LET k123[] = {1, 2, 3};
        LET k321[] = {3, 2, 1};

        TENSOR input_tensor;

        input_tensor.add_scal_prod(make_key(k11, 2), 1.0);
        input_tensor.add_scal_prod(make_key(k12, 2), 2.0);
        input_tensor.add_scal_prod(make_key(k21, 2), 3.0);
        input_tensor.add_scal_prod(make_key(k22, 2), 4.0);

        input_tensor.add_scal_prod(make_key(k123, 3), 5.0);

        TENSOR expected;

        expected.add_scal_prod(make_key(k11, 2), 1.0);
        expected.add_scal_prod(make_key(k12, 2), 3.0);
        expected.add_scal_prod(make_key(k21, 2), 2.0);
        expected.add_scal_prod(make_key(k22, 2), 4.0);

        expected.add_scal_prod(make_key(k321, 3), -5.0);

        TENSOR result;

        result = antipode(input_tensor);

        CHECK_EQUAL(expected, result);
    }

}// SUITE antipode
