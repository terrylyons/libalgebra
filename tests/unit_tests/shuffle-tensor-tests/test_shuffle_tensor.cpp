#include <iostream>

#include <UnitTest++/UnitTest++.h>

#include <libalgebra/libalgebra.h>


#include <libalgebra/alg_types.h>
#include <libalgebra/multiplication_helpers.h>



#include <libalgebra/tensor.h>

using alg::LET;
using alg::DEG;
using alg::shuffle_tensor;
using alg::shuffle_tensor_multiplication;


SUITE(shuffle_tensor)
{ 

    struct Fixture : public alg_types<5, 5, Rational>
    {
        typedef alg_types<5, 5, Rational> ALG_TYPES;
        typedef typename ALG_TYPES::TENSOR TENSOR;
        typedef typename TENSOR::BASIS TBASIS;
        typedef typename TENSOR::KEY KEY;

        using SHUFFLE_TENSOR_MULTIPLICATION = alg::shuffle_tensor_multiplication<typename ALG_TYPES::COEFF>;
        typedef typename ALG_TYPES::SHUFFLE_TENSOR SHUFFLE_TENSOR;

        // const TENSOR tunit;
        // const TENSOR tzero;

        const SHUFFLE_TENSOR tunit;
        const SHUFFLE_TENSOR tzero;

        Fixture() : tunit(KEY()), tzero()
        {}

        KEY make_key(const LET* arg, const std::size_t N)
        {
            KEY k;
            for (std::size_t i=0; i<N; ++i) {
                k.push_back(arg[i]);
            }
            return k;
        }

    };

    template <typename Coeff, DEG Width, DEG Depth>
    struct pairing{   
        
        using scalar_t = typename Coeff::S;
        using free_tensor_t = alg::free_tensor<Coeff, Width, Depth>;
        using shuffle_tensor_t = alg::shuffle_tensor<Coeff, Width, Depth>;

        scalar_t operator()(const shuffle_tensor_t& functional, const free_tensor_t* vector) const
        {
           scalar_t result {0};
           for (auto cit = functional.begin(); cit != functional.end(); ++cit) {
                result += cit->value() * (*vector)[cit->key()];
           }
           return result;
        }

    };

// #define ADD_KEY(N, ...) \
//     {                     \
//         LET tmp[N] = {__VA_ARGS__};  \
//         expected += TENSOR(make_key(tmp, N));\
//     }

    // test: {1()} * {} == {}

    TEST_FIXTURE(Fixture, test_star_tunit_zero) {


        // test: {1()} * {} == {}

        SHUFFLE_TENSOR lhs(tunit); // {1()}
        SHUFFLE_TENSOR rhs; // {}
        SHUFFLE_TENSOR expected; // {}

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        result = lhs * rhs;
            
                CHECK_EQUAL(expected, result);

    } // TEST test_star_tunit_zero

    TEST_FIXTURE(Fixture, test_multiply_and_add_tunit_zero_op) {

        // test: {1()} * {} == {}

        SHUFFLE_TENSOR lhs(tunit); // {1()}
        SHUFFLE_TENSOR rhs; // {}
        SHUFFLE_TENSOR expected; // {}

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough; 
        my_shuffle_tensor_product.multiply_and_add<SHUFFLE_TENSOR, alg::mult::scalar_passthrough>(result, lhs, rhs, my_scalar_passthrough);
            
                CHECK_EQUAL(expected, result);

    } // TEST test_multiply_and_add_tunit_zero_op

    TEST_FIXTURE(Fixture, test_multiply_and_add_tunit_zero_op_max_depth) {


        // test: {1()} * {} == {}

        SHUFFLE_TENSOR lhs(tunit); // {1()}
        SHUFFLE_TENSOR rhs; // {}
        SHUFFLE_TENSOR expected; // {}

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough; 
        my_shuffle_tensor_product.multiply_and_add<SHUFFLE_TENSOR, alg::mult::scalar_passthrough>(result, lhs, rhs, my_scalar_passthrough, 5);
            
                CHECK_EQUAL(expected, result);

    } // TEST test_multiply_and_add_tunit_zero_op_max_depth

    TEST_FIXTURE(Fixture, test_multiply_tunit_zero_op) {


        // test: {1()} * {} == {}

        SHUFFLE_TENSOR lhs(tunit); // {1()}
        SHUFFLE_TENSOR rhs; // {}
        SHUFFLE_TENSOR expected; // {}

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough; 
        result = my_shuffle_tensor_product.multiply<SHUFFLE_TENSOR, alg::mult::scalar_passthrough>(lhs, rhs, my_scalar_passthrough);
            
                CHECK_EQUAL(expected, result);

    } // TEST test_multiply_tunit_zero_op
    
    TEST_FIXTURE(Fixture, test_multiply_tunit_zero_op_max_depth) {


        // test: {1()} * {} == {}

        SHUFFLE_TENSOR lhs(tunit); // {1()}
        SHUFFLE_TENSOR rhs; // {}
        SHUFFLE_TENSOR expected; // {}

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough; 
        result = my_shuffle_tensor_product.multiply<SHUFFLE_TENSOR, alg::mult::scalar_passthrough>(lhs, rhs, my_scalar_passthrough, 5);
            
                CHECK_EQUAL(expected, result);

    } // TEST test_multiply_tunit_zero_op_max_depth

    TEST_FIXTURE(Fixture, test_multiply_inplace_tunit_zero_op) {


        // test: {1()} * {} == {}

        SHUFFLE_TENSOR lhs(tunit); // {1()}
        SHUFFLE_TENSOR rhs; // {}
        SHUFFLE_TENSOR expected; // {}

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough; 
        my_shuffle_tensor_product.multiply_inplace<SHUFFLE_TENSOR, alg::mult::scalar_passthrough>(lhs, rhs, my_scalar_passthrough);
            
                CHECK_EQUAL(expected, lhs);

    } // TEST test_multiply_inplace_tunit_zero_op
    
    TEST_FIXTURE(Fixture, test_multiply_inplace_tunit_zero_op_max_depth) {

        // test: {1()} * {} == {}

        SHUFFLE_TENSOR lhs(tunit); // {1()}
        SHUFFLE_TENSOR rhs; // {}
        SHUFFLE_TENSOR expected; // {}

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        my_shuffle_tensor_product.multiply_inplace<SHUFFLE_TENSOR, alg::mult::scalar_passthrough>(lhs, rhs, my_scalar_passthrough, 5); // TODO: Fix, does not compile
            
                CHECK_EQUAL(expected, lhs);

    } // TEST test_multiply_inplace_tunit_zero_op_max_depth

    // test: {} * {1()} == {}

    TEST_FIXTURE(Fixture, test_star_zero_tunit) {


        // test: {} * {1()} == {}

        SHUFFLE_TENSOR lhs; // {}
        SHUFFLE_TENSOR rhs(tunit); // {1()}
        SHUFFLE_TENSOR expected; // {}

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        result = lhs * rhs;

                CHECK_EQUAL(expected, result);

    } // TEST test_star_zero_tunit

    TEST_FIXTURE(Fixture, test_multiply_zero_tunit_op) {


        // test: {} * {1()} == {}

        SHUFFLE_TENSOR lhs; // {}
        SHUFFLE_TENSOR rhs(tunit); // {1()}
        SHUFFLE_TENSOR expected; // {}

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough; 
        result = my_shuffle_tensor_product.multiply<SHUFFLE_TENSOR, alg::mult::scalar_passthrough>(lhs, rhs, my_scalar_passthrough);

                CHECK_EQUAL(expected, result);

    } // TEST test_multiply_zero_tunit_op

    TEST_FIXTURE(Fixture, test_multiply_zero_tunit_op_max_depth) {


        // test: {} * {1()} == {}

        SHUFFLE_TENSOR lhs; // {}
        SHUFFLE_TENSOR rhs(tunit); // {1()}
        SHUFFLE_TENSOR expected; // {}

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough; 
        result = my_shuffle_tensor_product.multiply<SHUFFLE_TENSOR, alg::mult::scalar_passthrough>(lhs, rhs, my_scalar_passthrough, 5);

                CHECK_EQUAL(expected, result);

    } // TEST test_multiply_zero_tunit_op_max_depth

    // test: {1()} * {1()} == {1()}

    TEST_FIXTURE(Fixture, test_star_tunit_tunit) {


        // test: {1()} * {1()} == {1()}

        SHUFFLE_TENSOR lhs(tunit); // {1()}
        SHUFFLE_TENSOR rhs(tunit); // {1()}
        SHUFFLE_TENSOR expected(tunit); // {1()}

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        result = lhs * rhs;
            
                CHECK_EQUAL(expected, result);

    } // TEST test_star_tunit_tunit

    TEST_FIXTURE(Fixture, test_multiply_tunit_tunit_op) {


        // test: {1()} * {1()} == {1()}

        SHUFFLE_TENSOR lhs(tunit); // {1()}
        SHUFFLE_TENSOR rhs(tunit); // {1()}
        SHUFFLE_TENSOR expected(tunit); // {1()}

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough; 
        result = my_shuffle_tensor_product.multiply<SHUFFLE_TENSOR, alg::mult::scalar_passthrough>(lhs, rhs, my_scalar_passthrough);
            
                CHECK_EQUAL(expected, result);

    } // TEST test_multiply_tunit_tunit_op

    TEST_FIXTURE(Fixture, test_multiply_tunit_tunit_op_max_depth) {


    // test: {1()} * {1()} == {1()}

        SHUFFLE_TENSOR lhs(tunit); // {1()}
        SHUFFLE_TENSOR rhs(tunit); // {1()}
        SHUFFLE_TENSOR expected(tunit); // {1()}

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough; 
        result = my_shuffle_tensor_product.multiply<SHUFFLE_TENSOR, alg::mult::scalar_passthrough>(lhs, rhs, my_scalar_passthrough, 5);
            
                CHECK_EQUAL(expected, result);

    } // TEST test_multiply_tunit_tunit_op_max_depth

    // test: {1()} * {1(1)} ==  {1(1)}

    TEST_FIXTURE(Fixture, test_star_deg_1_tunit) {

        // test: {1()} * {1(1)} ==  {1(1)}

        SHUFFLE_TENSOR lhs(tunit); // {1()}

        LET k1[] = {1};

        SHUFFLE_TENSOR rhs(make_key(k1, 1)); // {1(1)}

        SHUFFLE_TENSOR expected(make_key(k1, 1)); // {1(1)}

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;

        result = lhs * rhs;

                CHECK_EQUAL(expected, result);

    } // TEST test_star_unidim_deg_1_tunit

    TEST_FIXTURE(Fixture, test_multiply_deg_1_tunit_op) {


        // test: {1()} * {1(1)} ==  {1(1)}

        SHUFFLE_TENSOR lhs(tunit); // {1()}

        LET k1[] = {1};

        SHUFFLE_TENSOR rhs(make_key(k1, 1)); // {1(1)}

        SHUFFLE_TENSOR expected(make_key(k1, 1)); // {1(1)}

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough; 
        result = my_shuffle_tensor_product.multiply<SHUFFLE_TENSOR, alg::mult::scalar_passthrough>(lhs, rhs, my_scalar_passthrough);

                CHECK_EQUAL(expected, result);

    } // TEST test_multiply_deg_1_tunit_op

    TEST_FIXTURE(Fixture, test_multiply_deg_1_tunit_op_max_depth) {

        // test: {1()} * {1(1)} ==  {1(1)}

        SHUFFLE_TENSOR lhs(tunit); // {1()}

        LET k1[] = {1};

        SHUFFLE_TENSOR rhs(make_key(k1, 1)); // {1(1)}

        SHUFFLE_TENSOR expected(make_key(k1, 1)); // {1(1)}

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough; 
        result = my_shuffle_tensor_product.multiply<SHUFFLE_TENSOR, alg::mult::scalar_passthrough>(lhs, rhs, my_scalar_passthrough, 5);

                CHECK_EQUAL(expected, result);

    } // TEST test_multiply_deg_1_tunit_op_max_depth

    // test: {1(1)} * {1(2)} == {1(1,2) 1(2,1)} 

    TEST_FIXTURE(Fixture, test_star_deg_1_deg_1) {

        // test: {1(1)} * {1(2)} == {1(1,2) 1(2,1)} 

        LET k1[] = {1};
        LET k2[] = {2};

        SHUFFLE_TENSOR lhs(make_key(k1, 1)); 
        SHUFFLE_TENSOR rhs(make_key(k2, 1));

        SHUFFLE_TENSOR expected;

        LET k12[] = {1, 2};
        LET k21[] = {2, 1};

        expected.add_scal_prod(make_key(k12,2), 1.0);
        expected.add_scal_prod(make_key(k21,2), 1.0);

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;

        result = lhs*rhs;

                CHECK_EQUAL(expected, result);

    } // TEST test_shuffle_product_unidim_deg_1_1

    TEST_FIXTURE(Fixture, test_multiply_deg_1_deg_1_op) {

        // test: {1(1)} * {1(2)} == {1(1,2) 1(2,1)}

        LET k1[] = {1};
        LET k2[] = {2};

        SHUFFLE_TENSOR lhs(make_key(k1, 1)); 
        SHUFFLE_TENSOR rhs(make_key(k2, 1));

        SHUFFLE_TENSOR expected;

        LET k12[] = {1, 2};
        LET k21[] = {2, 1};

        expected.add_scal_prod(make_key(k12,2), 1.0);
        expected.add_scal_prod(make_key(k21,2), 1.0);

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough; 
        result = my_shuffle_tensor_product.multiply<SHUFFLE_TENSOR, alg::mult::scalar_passthrough>(lhs, rhs, my_scalar_passthrough);

                CHECK_EQUAL(expected, result);

    } // TEST test_multiply_deg_1_deg_1_op

    TEST_FIXTURE(Fixture, test_multiply_deg_1_deg_1_op_max_depth) {

        // test: {1(1)} * {1(2)} == {1(1,2) 1(2,1)}

        LET k1[] = {1};
        LET k2[] = {2};

        SHUFFLE_TENSOR lhs(make_key(k1, 1)); 
        SHUFFLE_TENSOR rhs(make_key(k2, 1));

        SHUFFLE_TENSOR expected;

        LET k12[] = {1, 2};
        LET k21[] = {2, 1};

        expected.add_scal_prod(make_key(k12,2), 1.0);
        expected.add_scal_prod(make_key(k21,2), 1.0);

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough; 
        result = my_shuffle_tensor_product.multiply<SHUFFLE_TENSOR, alg::mult::scalar_passthrough>(lhs, rhs, my_scalar_passthrough, 5);

                CHECK_EQUAL(expected, result);

    } // TEST test_multiply_deg_1_deg_1_op_max_depth

    // test: {1(1)} * {1(2,3)} == {1(1,2,3) 1(2,1,3) 1(2,3,1)}

    TEST_FIXTURE(Fixture, test_star_deg_1_deg_2) {

        // test: {1(1)} * {1(2,3)} == {1(1,2,3) 1(2,1,3) 1(2,3,1)}

        LET k1[] = {1};
        LET k2[] = {2, 3};

        SHUFFLE_TENSOR lhs(make_key(k1, 1)); // {1()1}
        SHUFFLE_TENSOR rhs(make_key(k2, 2)); // {1(2,3)}

        SHUFFLE_TENSOR expected;

        LET k123[] = {1, 2, 3};
        LET k213[] = {2, 1, 3};
        LET k231[] = {2, 3, 1};

        expected.add_scal_prod(make_key(k123,3), 1.0);
        expected.add_scal_prod(make_key(k213,3), 1.0);
        expected.add_scal_prod(make_key(k231,3), 1.0);

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;

        result = lhs * rhs;

                CHECK_EQUAL(expected, result);

    } // TEST test_star_unidim_deg_1_deg_2

    TEST_FIXTURE(Fixture, test_multiply_deg_1_deg_2_op) {

        // test: {1(1)} * {1(2,3)} == {1(1,2,3) 1(2,1,3) 1(2,3,1)}

        LET k1[] = {1};
        LET k2[] = {2, 3};

        SHUFFLE_TENSOR lhs(make_key(k1, 1)); // {1()1}
        SHUFFLE_TENSOR rhs(make_key(k2, 2)); // {1(2,3)}

        SHUFFLE_TENSOR expected;

        LET k123[] = {1, 2, 3};
        LET k213[] = {2, 1, 3};
        LET k231[] = {2, 3, 1};

        expected.add_scal_prod(make_key(k123,3), 1.0);
        expected.add_scal_prod(make_key(k213,3), 1.0);
        expected.add_scal_prod(make_key(k231,3), 1.0);

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough; 
        result = my_shuffle_tensor_product.multiply<SHUFFLE_TENSOR, alg::mult::scalar_passthrough>(lhs, rhs, my_scalar_passthrough);

                CHECK_EQUAL(expected, result);

    } // TEST test_multiply_deg_1_deg_2_op

    TEST_FIXTURE(Fixture, test_multiply_deg_1_deg_2_op_max_depth) {

        // test: {1(1)} * {1(2,3)} == {1(1,2,3) 1(2,1,3) 1(2,3,1)}

        LET k1[] = {1};
        LET k2[] = {2, 3};

        SHUFFLE_TENSOR lhs(make_key(k1, 1)); // {1()1}
        SHUFFLE_TENSOR rhs(make_key(k2, 2)); // {1(2,3)}

        SHUFFLE_TENSOR expected;

        LET k123[] = {1, 2, 3};
        LET k213[] = {2, 1, 3};
        LET k231[] = {2, 3, 1};

        expected.add_scal_prod(make_key(k123,3), 1.0);
        expected.add_scal_prod(make_key(k213,3), 1.0);
        expected.add_scal_prod(make_key(k231,3), 1.0);

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough; 
        result = my_shuffle_tensor_product.multiply<SHUFFLE_TENSOR, alg::mult::scalar_passthrough>(lhs, rhs, my_scalar_passthrough, 5);

                CHECK_EQUAL(expected, result);

    } // TEST test_multiply_deg_1_deg_2_op_max_depth

    // test: {1(1,2)} * {1(3)} == {1(1,2,3) 1(1,3,2) 1(3,1,2)}

    TEST_FIXTURE(Fixture, test_star_deg_2_deg_1) {

        // test: {1(1,2)} * {1(3)} == {1(1,2,3) 1(1,3,2) 1(3,1,2)}

        LET k1[] = {1, 2};
        LET k2[] = {3};

        SHUFFLE_TENSOR lhs(make_key(k1, 2)); // {1(1,2)}
        SHUFFLE_TENSOR rhs(make_key(k2, 1)); // {1(3)}

        SHUFFLE_TENSOR expected;

        LET k123[] = {1, 2, 3};
        LET k132[] = {1, 3, 2};
        LET k312[] = {3, 1, 2};

        expected.add_scal_prod(make_key(k123,3), 1.0); // 1(1,2,3)
        expected.add_scal_prod(make_key(k132,3), 1.0); // 1(1,3,2)
        expected.add_scal_prod(make_key(k312,3), 1.0); // 1(3,1,2)      // {1(1,2,3) 1(1,3,2) 1(3,1,2)}

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;

        result = lhs * rhs;

                CHECK_EQUAL(expected, result);
            
    } // TEST test_product_unidim_deg_2_deg_1

    TEST_FIXTURE(Fixture, test_multiply_deg_2_deg_1_op) {

        // test: {1(1,2)} * {1(3)} == {1(1,2,3) 1(1,3,2) 1(3,1,2)}

        LET k1[] = {1, 2};
        LET k2[] = {3};

        SHUFFLE_TENSOR lhs(make_key(k1, 2)); // {1(1,2)}
        SHUFFLE_TENSOR rhs(make_key(k2, 1)); // {1(3)}

        SHUFFLE_TENSOR expected;

        LET k123[] = {1, 2, 3};
        LET k132[] = {1, 3, 2};
        LET k312[] = {3, 1, 2};

        expected.add_scal_prod(make_key(k123,3), 1.0); // 1(1,2,3)
        expected.add_scal_prod(make_key(k132,3), 1.0); // 1(1,3,2)
        expected.add_scal_prod(make_key(k312,3), 1.0); // 1(3,1,2)      // {1(1,2,3) 1(1,3,2) 1(3,1,2)}

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough; 
        result = my_shuffle_tensor_product.multiply<SHUFFLE_TENSOR, alg::mult::scalar_passthrough>(lhs, rhs, my_scalar_passthrough);

                CHECK_EQUAL(expected, result);
            
    } // TEST test_multiply_deg_2_deg_1_op

    TEST_FIXTURE(Fixture, test_multiply_deg_2_deg_1_op_max_depth) {

        // test: {1(1,2)} * {1(3)} == {1(1,2,3) 1(1,3,2) 1(3,1,2)}

        LET k1[] = {1, 2};
        LET k2[] = {3};

        SHUFFLE_TENSOR lhs(make_key(k1, 2)); // {1(1,2)}
        SHUFFLE_TENSOR rhs(make_key(k2, 1)); // {1(3)}

        SHUFFLE_TENSOR expected;

        LET k123[] = {1, 2, 3};
        LET k132[] = {1, 3, 2};
        LET k312[] = {3, 1, 2};

        expected.add_scal_prod(make_key(k123,3), 1.0); // 1(1,2,3)
        expected.add_scal_prod(make_key(k132,3), 1.0); // 1(1,3,2)
        expected.add_scal_prod(make_key(k312,3), 1.0); // 1(3,1,2)      // {1(1,2,3) 1(1,3,2) 1(3,1,2)}

        std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough; 
        result = my_shuffle_tensor_product.multiply<SHUFFLE_TENSOR, alg::mult::scalar_passthrough>(lhs, rhs, my_scalar_passthrough, 5);

                CHECK_EQUAL(expected, result);
            
    } // TEST test_multiply_deg_2_deg_1_op_max_depth

} // SUITE shuffle_tensor
