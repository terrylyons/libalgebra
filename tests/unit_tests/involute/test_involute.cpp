#include <UnitTest++/UnitTest++.h>

#include <libalgebra/libalgebra.h>

#include <libalgebra/alg_types.h>
#include <libalgebra/multiplication_helpers.h>

#include <libalgebra/tensor.h>

#include "../../common/random_vector_generator.h"
#include "../../common/rng.h"

using alg::DEG;
using alg::LET;

SUITE(involute)
{
    struct Fixture : public alg_types<5, 5, Rational> {
        typedef alg_types<5, 5, Rational> ALG_TYPES;
        typedef typename ALG_TYPES::TENSOR TENSOR;
        typedef typename TENSOR::BASIS TBASIS;
        typedef typename TENSOR::KEY KEY;

         const TENSOR tunit;
         const TENSOR tzero;

         Fixture() : tunit(KEY()), tzero()
         {}

//         TENSOR reference_tensor;
//         alg::operators::left_multiplication_operator<TENSOR> op;
//
//         using rat_dist = la_testing::uniform_rational_distribution<S>;
//         using rvg_t = la_testing::random_vector_generator<TENSOR, rat_dist>;
//
//         std::mt19937 rng;
//         rvg_t rvg;
//
//         Fixture() : reference_tensor(typename TENSOR::KEY(alg::LET(1))),
//                     op(TENSOR(reference_tensor)), rng(std::random_device()()), rvg(-1, 1),
//                     tunit(KEY()), tzero()
//         {}

//         KEY make_key(const LET* arg, const std::size_t N)
//         {
//             KEY k;
//             for (std::size_t i = 0; i < N; ++i) {
//                 k.push_back(arg[i]);
//             }
//             return k;
//         }
    };

//    TEST_FIXTURE(Fixture, UnitTest1)
//    {
//
//        CHECK_EQUAL(0, 0);
//
//    }

//    TEST_FIXTURE(Fixture, UnitTest2) {
//
//        TENSOR my_random_tensor = rvg(rng);
//
////        std::cout << "my_random_tensor=" << my_random_tensor << std::endl;
//
//        CHECK_EQUAL(0, 0);
//    }

    TEST_FIXTURE(Fixture, UnitTest) {

        TENSOR my_tensor(tunit);

//        alg::vectors::dtl::data_access_base<> my_data_access_base;

        CHECK_EQUAL(0, 1);
    }

}// SUITE involute