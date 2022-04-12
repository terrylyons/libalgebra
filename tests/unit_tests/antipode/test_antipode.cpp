#include <UnitTest++/UnitTest++.h>

#include <libalgebra/libalgebra.h>

#include <iostream>
#include <vector>

#include <libalgebra/alg_types.h>
#include <libalgebra/multiplication_helpers.h>

#include <libalgebra/operators/operators.h>

#include "../../common/random_vector_generator.h"
#include "../../common/rng.h"

#include <libalgebra/tensor.h>

using alg::LET;

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
    typedef alg::lie<Coeff, Width, Depth, alg::vectors::dense_vector> LIE;

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

struct RandomFixture : alg_types<5, 5, Rational> {
    TENSOR reference_tensor;
    alg::operators::left_multiplication_operator<TENSOR> op;

    using rat_dist = la_testing::uniform_rational_distribution<S>;
    using rvg_t = la_testing::random_vector_generator<TENSOR, rat_dist>;

    std::mt19937 rng;
    rvg_t rvg;

    const TENSOR tunit;

    typedef typename TENSOR::KEY KEY;
    KEY kunit;

    RandomFixture() : reference_tensor(typename TENSOR::KEY(alg::LET(1))),
                      op(TENSOR(reference_tensor)), rng(std::random_device()()), rvg(-1, 1)
    {}
};

SUITE(Antipode)
{
    typedef alg::coefficients::coefficient_field<float> float_field;
    typedef DenseFixture<float_field, 4, 4> dense_fixture;

    // test: antipode(zero) == zero

    TEST_FIXTURE(dense_fixture, DenseAntipodeZero)
    {
        // check: { } --> { }

        TENSOR result = antipode(tzero);

        CHECK_EQUAL(tzero, result);
    }

    // test: antipode(identity) == identity

    TEST_FIXTURE(dense_fixture, DenseAntipodeIdentiy)
    {
        // check: { 1{} } --> { 1{} }

        TENSOR result = antipode(tunit);

        CHECK_EQUAL(tunit, result);
    }

    // test: antipode(length 1 word) == - length 1 word

    TEST_FIXTURE(dense_fixture, DenseAntipodeOneLetter)
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

    TEST_FIXTURE(dense_fixture, DenseAntipodeOneEvenWord)
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

    TEST_FIXTURE(dense_fixture, DenseAntipodeOneOddWord)
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

    TEST_FIXTURE(dense_fixture, DenseAntipodeMultipleWord)
    {
        // check: { 1.0{11} 2.0{12} 3.0{21} 4.0{22} } --> { 1.0{11} 3.0{12} 2.0{21} 4.0{22} }

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

        std::cout << "input_tensor=" << input_tensor << std::endl;

        std::cout << "expected=" << expected << std::endl;

        std::cout << "result=" << result << std::endl;

        CHECK_EQUAL(expected, result);
    }

    // tests using random tensors: antipode(a) * a == Identity

    TEST_FIXTURE(RandomFixture, DenseAntipodeGroupLike)
    {
        LIE incr;
        incr.add_scal_prod(typename LIE::KEY(1), typename LIE::SCALAR(1));
        incr.add_scal_prod(typename LIE::KEY(2), typename LIE::SCALAR(2));
        incr.add_scal_prod(typename LIE::KEY(3), typename LIE::SCALAR(3));
        incr.add_scal_prod(typename LIE::KEY(4), typename LIE::SCALAR(4));

        alg::maps<dense_fixtureDenseAntipodeIdentiyHelper::coeffs, dense_fixtureDenseAntipodeIdentiyHelper::width, dense_fixtureDenseAntipodeIdentiyHelper::depth, TENSOR, LIE> maps;

        auto input_tensor = exp(maps.l2t(incr));

        auto result = antipode(input_tensor) * input_tensor;

        double tolerance = 0.01;// TODO: change this!

        for (auto cit = result.begin(); cit != result.end(); ++cit) {
            if (cit->key() == kunit) {
                // check 1
                CHECK_CLOSE(1.0, cit->value(), tolerance);
            }
            else {
                // check 0
                CHECK_CLOSE(0.0, cit->value(), tolerance);
            }
        }
    }

    typedef SparseFixture<float_field, 4, 4> sparse_fixture;

    // test: key/index look-ups for multiple word tensor

    TEST_FIXTURE(sparse_fixture, SparseAntipodeMultipleWord)
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

        //            std::cout << "input_tensor=" << input_tensor << std::endl;
        //
        //            std::cout << "expected=" << expected << std::endl;
        //
        //            std::cout << "result=" << result << std::endl;

        CHECK_EQUAL(expected, result);
    }

}// SUITE antipode