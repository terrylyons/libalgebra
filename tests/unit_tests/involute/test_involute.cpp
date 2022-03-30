#include <UnitTest++/UnitTest++.h>


#include <libalgebra/libalgebra.h>

#include <vector>

#include <libalgebra/alg_types.h>
#include <libalgebra/multiplication_helpers.h>

#include <libalgebra/tensor.h>

using alg::LET;

template <typename Coeff, unsigned Width, unsigned Depth>
struct SparseFixture
{

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


template <typename Coeff, unsigned Width, unsigned Depth>
struct DenseFixture
{

    typedef typename Coeff::S S;
    typedef typename Coeff::Q Q;

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


SUITE(involute)
{
        typedef alg::coefficients::coefficient_field<float> float_field;
        typedef DenseFixture<float_field, 4, 4> dense_fixture;

        // involute({ 1.0{11} 2.0{12} 3.0{21} 4.0{22} }) == { 1.0{11} 3.0{12} 2.0{21} 4.0{22} }

        TEST_FIXTURE(dense_fixture, dense_unit_test)
        {

            LET k11[] = {1, 1};
            LET k12[] = {1, 2};
            LET k21[] = {2, 1};
            LET k22[] = {2, 2};

            TENSOR input_tensor;

            input_tensor.add_scal_prod(make_key(k11, 2), 1.0);
            input_tensor.add_scal_prod(make_key(k12, 2), 2.0);
            input_tensor.add_scal_prod(make_key(k21, 2), 3.0);
            input_tensor.add_scal_prod(make_key(k22, 2), 4.0);

            TENSOR expected;

            expected.add_scal_prod(make_key(k11, 2), 1.0);
            expected.add_scal_prod(make_key(k12, 2), 3.0);
            expected.add_scal_prod(make_key(k21, 2), 2.0);
            expected.add_scal_prod(make_key(k22, 2), 4.0);

            TENSOR result;

            result = involute(input_tensor);

            std::cout << "input_tensor=" << input_tensor << std::endl;

            std::cout << "expected=" << expected << std::endl;

            std::cout << "result=" << result << std::endl;

            CHECK_EQUAL(expected, result);
        }

        typedef SparseFixture<float_field, 4, 4> sparse_fixture;

        // involute({ 1.0{11} 2.0{12} 3.0{21} 4.0{22} 5.0{123}}) == { 1.0{11} 3.0{12} 2.0{21} 4.0{22} -5.0{321}}

        TEST_FIXTURE(sparse_fixture, sparse_unit_test)
        {

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

            result = involute(input_tensor);

//            std::cout << "input_tensor=" << input_tensor << std::endl;
//
//            std::cout << "expected=" << expected << std::endl;
//
//            std::cout << "result=" << result << std::endl;

            CHECK_EQUAL(expected, result);
        }


}// SUITE involute