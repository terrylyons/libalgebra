#include <UnitTest++/UnitTest++.h>


#include <libalgebra/libalgebra.h>

#include <vector>

#include <libalgebra/alg_types.h>
#include <libalgebra/multiplication_helpers.h>

#include <libalgebra/tensor.h>

using alg::LET;


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
        typedef DenseFixture<float_field, 5, 5> dense_fixture;

        // involute({12(1,2) 21(2,1)}) == {21(1,2) 12(2,1)}

        TEST_FIXTURE(dense_fixture, dense_unit_test)
        {

            LET k12[] = {1, 2};
            LET k21[] = {2, 1};

            TENSOR input_tensor;

            input_tensor.add_scal_prod(make_key(k12, 2), 12.0);
            input_tensor.add_scal_prod(make_key(k21, 2), 21.0);

            TENSOR expected;

            expected.add_scal_prod(make_key(k12, 2), 21.0);
            expected.add_scal_prod(make_key(k21, 2), 12.0);

            TENSOR result;

            result = involute(input_tensor);

            std::cout << "input_tensor=" << input_tensor << std::endl;

            std::cout << "expected=" << expected << std::endl;

            std::cout << "result=" << result << std::endl;

            CHECK_EQUAL(expected, result);
        }


}// SUITE involute