#include <iostream>

#include <UnitTest++.h>

#include <libalgebra/libalgebra.h>

#include <libalgebra/alg_types.h>
#include <libalgebra/multiplication_helpers.h>

#include <libalgebra/tensor.h>

using alg::DEG;
using alg::LET;
using alg::shuffle_tensor;
using alg::shuffle_tensor_multiplication;

SUITE(shuffle_tensor)
{

    struct Fixture : public alg_types<5, 5, Rational> {
        typedef alg_types<5, 5, Rational> ALG_TYPES;
        typedef typename ALG_TYPES::TENSOR TENSOR;
        typedef typename TENSOR::BASIS TBASIS;
        typedef typename TENSOR::KEY KEY;

        using SHUFFLE_TENSOR_MULTIPLICATION = alg::shuffle_tensor_multiplication<5, 5>;
        typedef typename ALG_TYPES::SHUFFLE_TENSOR SHUFFLE_TENSOR;
        using mtraits = alg::dtl::multiplication_traits<SHUFFLE_TENSOR_MULTIPLICATION>;
        // const TENSOR tunit;
        // const TENSOR tzero;

        const SHUFFLE_TENSOR tunit;
        const SHUFFLE_TENSOR tzero;

        Fixture() : tunit(KEY()), tzero()
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

    template<typename Coeff, DEG Width, DEG Depth>
    struct pairing {

        using scalar_t = typename Coeff::S;
        using free_tensor_t = alg::free_tensor<Coeff, Width, Depth>;
        using shuffle_tensor_t = alg::shuffle_tensor<Coeff, Width, Depth>;

        scalar_t operator()(const shuffle_tensor_t& functional, const free_tensor_t& vector) const
        {
            scalar_t result{0};
            for (auto cit = functional.begin(); cit != functional.end(); ++cit) {
                result += cit->value() * vector[cit->key()];
            }
            return result;
        }
    };

    // test: {1()} * {} == {}

    TEST_FIXTURE(Fixture, test_star_tunit_zero)
    {

        // test: {1()} * {} == {}

        SHUFFLE_TENSOR lhs(tunit);// {1()}
        SHUFFLE_TENSOR rhs;       // {}
        SHUFFLE_TENSOR expected;  // {}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        result = lhs * rhs;

        CHECK_EQUAL(expected, result);

    }// TEST test_star_tunit_zero

    TEST_FIXTURE(Fixture, test_multiply_and_add_tunit_zero_op)
    {

        // test: {1()} * {} == {}

        SHUFFLE_TENSOR lhs(tunit);// {1()}
        SHUFFLE_TENSOR rhs;       // {}
        SHUFFLE_TENSOR expected;  // {}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;

        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_and_add_tunit_zero_op

    TEST_FIXTURE(Fixture, test_multiply_and_add_tunit_zero_op_max_depth)
    {

        // test: {1()} * {} == {}

        SHUFFLE_TENSOR lhs(tunit);// {1()}
        SHUFFLE_TENSOR rhs;       // {}
        SHUFFLE_TENSOR expected;  // {}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;

        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough, 5);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_and_add_tunit_zero_op_max_depth

    TEST_FIXTURE(Fixture, test_multiply_tunit_zero_op)
    {

        // test: {1()} * {} == {}

        SHUFFLE_TENSOR lhs(tunit);// {1()}
        SHUFFLE_TENSOR rhs;       // {}
        SHUFFLE_TENSOR expected;  // {}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_tunit_zero_op

    TEST_FIXTURE(Fixture, test_multiply_tunit_zero_op_max_depth)
    {

        // test: {1()} * {} == {}

        SHUFFLE_TENSOR lhs(tunit);// {1()}
        SHUFFLE_TENSOR rhs;       // {}
        SHUFFLE_TENSOR expected;  // {}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough, 5);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_tunit_zero_op_max_depth

    TEST_FIXTURE(Fixture, test_multiply_inplace_tunit_zero_op)
    {

        // test: {1()} * {} == {}

        SHUFFLE_TENSOR lhs(tunit);// {1()}
        SHUFFLE_TENSOR rhs;       // {}
        SHUFFLE_TENSOR expected;  // {}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_inplace(my_shuffle_tensor_product, lhs, rhs, my_scalar_passthrough);

        CHECK_EQUAL(expected, lhs);

    }// TEST test_multiply_inplace_tunit_zero_op

    TEST_FIXTURE(Fixture, test_multiply_inplace_tunit_zero_op_max_depth)
    {

        // test: {1()} * {} == {}

        SHUFFLE_TENSOR lhs(tunit);// {1()}
        SHUFFLE_TENSOR rhs;       // {}
        SHUFFLE_TENSOR expected;  // {}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_inplace(my_shuffle_tensor_product, lhs, rhs, my_scalar_passthrough, 5);

        CHECK_EQUAL(expected, lhs);

    }// TEST test_multiply_inplace_tunit_zero_op_max_depth

    // test: {} * {1()} == {}

    TEST_FIXTURE(Fixture, test_star_zero_tunit)
    {

        // test: {} * {1()} == {}

        SHUFFLE_TENSOR lhs;       // {}
        SHUFFLE_TENSOR rhs(tunit);// {1()}
        SHUFFLE_TENSOR expected;  // {}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        result = lhs * rhs;

        CHECK_EQUAL(expected, result);

    }// TEST test_star_zero_tunit

    TEST_FIXTURE(Fixture, test_multiply_and_add_zero_tunit_op)
    {

        // test: {} * {1()} == {}

        SHUFFLE_TENSOR lhs;       // {}
        SHUFFLE_TENSOR rhs(tunit);// {1()}
        SHUFFLE_TENSOR expected;  // {}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_and_add_zero_tunit_op

    TEST_FIXTURE(Fixture, test_multiply_and_add_zero_tunit_op_max_depth)
    {

        // test: {} * {1()} == {}

        SHUFFLE_TENSOR lhs;       // {}
        SHUFFLE_TENSOR rhs(tunit);// {1()}
        SHUFFLE_TENSOR expected;  // {}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough, 5);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_and_add_zero_tunit_op_max_depth

    TEST_FIXTURE(Fixture, test_multiply_zero_tunit_op)
    {

        // test: {} * {1()} == {}

        SHUFFLE_TENSOR lhs;       // {}
        SHUFFLE_TENSOR rhs(tunit);// {1()}
        SHUFFLE_TENSOR expected;  // {}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_zero_tunit_op

    TEST_FIXTURE(Fixture, test_multiply_zero_tunit_op_max_depth)
    {

        // test: {} * {1()} == {}

        SHUFFLE_TENSOR lhs;       // {}
        SHUFFLE_TENSOR rhs(tunit);// {1()}
        SHUFFLE_TENSOR expected;  // {}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough, 5);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_zero_tunit_op_max_depth

    TEST_FIXTURE(Fixture, test_multiply_inplace_zero_tunit_op)
    {

        // test: {} * {1()} == {}

        SHUFFLE_TENSOR lhs;       // {}
        SHUFFLE_TENSOR rhs(tunit);// {1()}
        SHUFFLE_TENSOR expected;  // {}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_inplace(my_shuffle_tensor_product, lhs, rhs, my_scalar_passthrough);

        CHECK_EQUAL(expected, lhs);

    }// TEST test_multiply_inplace_zero_tunit_op

    TEST_FIXTURE(Fixture, test_multiply_inplace_zero_tunit_op_max_depth)
    {

        // test: {} * {1()} == {}

        SHUFFLE_TENSOR lhs;       // {}
        SHUFFLE_TENSOR rhs(tunit);// {1()}
        SHUFFLE_TENSOR expected;  // {}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_inplace(my_shuffle_tensor_product, lhs, rhs, my_scalar_passthrough, 5);

        CHECK_EQUAL(expected, lhs);

    }// TEST test_multiply_inplace_zero_tunit_op_max_depth

    // test: {1()} * {1()} == {1()}

    TEST_FIXTURE(Fixture, test_star_tunit_tunit)
    {

        // test: {1()} * {1()} == {1()}

        SHUFFLE_TENSOR lhs(tunit);     // {1()}
        SHUFFLE_TENSOR rhs(tunit);     // {1()}
        SHUFFLE_TENSOR expected(tunit);// {1()}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        result = lhs * rhs;

        CHECK_EQUAL(expected, result);

    }// TEST test_star_tunit_tunit

    TEST_FIXTURE(Fixture, test_multiply_and_add_tunit_tunit_op)
    {

        // test: {1()} * {1()} == {1()}

        SHUFFLE_TENSOR lhs(tunit);     // {1()}
        SHUFFLE_TENSOR rhs(tunit);     // {1()}
        SHUFFLE_TENSOR expected(tunit);// {1()}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_and_add_tunit_tunit_op

    TEST_FIXTURE(Fixture, test_multiply_and_add_tunit_tunit_op_max_depth)
    {

        // test: {1()} * {1()} == {1()}

        SHUFFLE_TENSOR lhs(tunit);     // {1()}
        SHUFFLE_TENSOR rhs(tunit);     // {1()}
        SHUFFLE_TENSOR expected(tunit);// {1()}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough, 5);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_and_add_tunit_tunit_op_max_depth

    TEST_FIXTURE(Fixture, test_multiply_tunit_tunit_op)
    {

        // test: {1()} * {1()} == {1()}

        SHUFFLE_TENSOR lhs(tunit);     // {1()}
        SHUFFLE_TENSOR rhs(tunit);     // {1()}
        SHUFFLE_TENSOR expected(tunit);// {1()}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_tunit_tunit_op

    TEST_FIXTURE(Fixture, test_multiply_tunit_tunit_op_max_depth)
    {

        // test: {1()} * {1()} == {1()}

        SHUFFLE_TENSOR lhs(tunit);     // {1()}
        SHUFFLE_TENSOR rhs(tunit);     // {1()}
        SHUFFLE_TENSOR expected(tunit);// {1()}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough, 5);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_tunit_tunit_op_max_depth

    TEST_FIXTURE(Fixture, test_multiply_inplace_tunit_tunit_op)
    {

        // test: {1()} * {1()} == {1()}

        SHUFFLE_TENSOR lhs(tunit);     // {1()}
        SHUFFLE_TENSOR rhs(tunit);     // {1()}
        SHUFFLE_TENSOR expected(tunit);// {1()}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_inplace(my_shuffle_tensor_product, lhs, rhs, my_scalar_passthrough);

        CHECK_EQUAL(expected, lhs);

    }// TEST test_multiply_inplace_tunit_tunit_op

    TEST_FIXTURE(Fixture, test_multiply_inplace_tunit_tunit_op_max_depth)
    {

        // test: {1()} * {1()} == {1()}

        SHUFFLE_TENSOR lhs(tunit);     // {1()}
        SHUFFLE_TENSOR rhs(tunit);     // {1()}
        SHUFFLE_TENSOR expected(tunit);// {1()}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_inplace(my_shuffle_tensor_product, lhs, rhs, my_scalar_passthrough, 5);

        CHECK_EQUAL(expected, lhs);

    }// TEST test_multiply_inplace_tunit_tunit_op_max_depth

    // test: {1()} * {1(1)} ==  {1(1)}

    TEST_FIXTURE(Fixture, test_star_deg_1_tunit)
    {

        // test: {1()} * {1(1)} ==  {1(1)}

        SHUFFLE_TENSOR lhs(tunit);// {1()}

        LET k1[] = {1};

        SHUFFLE_TENSOR rhs(make_key(k1, 1));// {1(1)}

        SHUFFLE_TENSOR expected(make_key(k1, 1));// {1(1)}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;

        result = lhs * rhs;

        CHECK_EQUAL(expected, result);

    }// TEST test_star_unidim_deg_1_tunit

    TEST_FIXTURE(Fixture, test_multiply_and_add_deg_1_tunit_op)
    {

        // test: {1()} * {1(1)} ==  {1(1)}

        SHUFFLE_TENSOR lhs(tunit);// {1()}

        LET k1[] = {1};

        SHUFFLE_TENSOR rhs(make_key(k1, 1));// {1(1)}

        SHUFFLE_TENSOR expected(make_key(k1, 1));// {1(1)}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_and_add_deg_1_tunit_op

    TEST_FIXTURE(Fixture, test_multiply_and_add_deg_1_tunit_op_max_depth)
    {

        // test: {1()} * {1(1)} ==  {1(1)}

        SHUFFLE_TENSOR lhs(tunit);// {1()}

        LET k1[] = {1};

        SHUFFLE_TENSOR rhs(make_key(k1, 1));// {1(1)}

        SHUFFLE_TENSOR expected(make_key(k1, 1));// {1(1)}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough, 5);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_and_add_deg_1_tunit_op_max_depth

    TEST_FIXTURE(Fixture, test_multiply_deg_1_tunit_op)
    {

        // test: {1()} * {1(1)} ==  {1(1)}

        SHUFFLE_TENSOR lhs(tunit);// {1()}

        LET k1[] = {1};

        SHUFFLE_TENSOR rhs(make_key(k1, 1));// {1(1)}

        SHUFFLE_TENSOR expected(make_key(k1, 1));// {1(1)}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_deg_1_tunit_op

    TEST_FIXTURE(Fixture, test_multiply_deg_1_tunit_op_max_depth)
    {

        // test: {1()} * {1(1)} ==  {1(1)}

        SHUFFLE_TENSOR lhs(tunit);// {1()}

        LET k1[] = {1};

        SHUFFLE_TENSOR rhs(make_key(k1, 1));// {1(1)}

        SHUFFLE_TENSOR expected(make_key(k1, 1));// {1(1)}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough, 5);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_deg_1_tunit_op_max_depth

    TEST_FIXTURE(Fixture, test_multiply_inplace_deg_1_tunit_op)
    {

        // test: {1()} * {1(1)} ==  {1(1)}

        SHUFFLE_TENSOR lhs(tunit);// {1()}

        LET k1[] = {1};

        SHUFFLE_TENSOR rhs(make_key(k1, 1));// {1(1)}

        SHUFFLE_TENSOR expected(make_key(k1, 1));// {1(1)}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_inplace(my_shuffle_tensor_product, lhs, rhs, my_scalar_passthrough);

        CHECK_EQUAL(expected, lhs);

    }// TEST test_multiply_inplace_deg_1_tunit_op

    TEST_FIXTURE(Fixture, test_multiply_inplace_deg_1_tunit_op_max_depth)
    {

        // test: {1()} * {1(1)} ==  {1(1)}

        SHUFFLE_TENSOR lhs(tunit);// {1()}

        LET k1[] = {1};

        SHUFFLE_TENSOR rhs(make_key(k1, 1));// {1(1)}

        SHUFFLE_TENSOR expected(make_key(k1, 1));// {1(1)}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_inplace(my_shuffle_tensor_product, lhs, rhs, my_scalar_passthrough, 5);

        CHECK_EQUAL(expected, lhs);

    }// TEST test_multiply_inplace_deg_1_tunit_op_max_depth

    // test: {1(1)} * {1(2)} == {1(1,2) 1(2,1)}

    TEST_FIXTURE(Fixture, test_star_deg_1_deg_1)
    {

        // test: {1(1)} * {1(2)} == {1(1,2) 1(2,1)}

        LET k1[] = {1};
        LET k2[] = {2};

        SHUFFLE_TENSOR lhs(make_key(k1, 1));
        SHUFFLE_TENSOR rhs(make_key(k2, 1));

        SHUFFLE_TENSOR expected;

        LET k12[] = {1, 2};
        LET k21[] = {2, 1};

        expected.add_scal_prod(make_key(k12, 2), 1.0);
        expected.add_scal_prod(make_key(k21, 2), 1.0);

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;

        result = lhs * rhs;

        CHECK_EQUAL(expected, result);

    }// TEST test_shuffle_product_unidim_deg_1_1

    TEST_FIXTURE(Fixture, test_multiply_and_add_deg_1_deg_1_op)
    {

        // test: {1(1)} * {1(2)} == {1(1,2) 1(2,1)}

        LET k1[] = {1};
        LET k2[] = {2};

        SHUFFLE_TENSOR lhs(make_key(k1, 1));
        SHUFFLE_TENSOR rhs(make_key(k2, 1));

        SHUFFLE_TENSOR expected;

        LET k12[] = {1, 2};
        LET k21[] = {2, 1};

        expected.add_scal_prod(make_key(k12, 2), 1.0);
        expected.add_scal_prod(make_key(k21, 2), 1.0);

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_and_add_deg_1_deg_1_op

    TEST_FIXTURE(Fixture, test_multiply_and_add_deg_1_deg_1_op_max_depth)
    {

        // test: {1(1)} * {1(2)} == {1(1,2) 1(2,1)}

        LET k1[] = {1};
        LET k2[] = {2};

        SHUFFLE_TENSOR lhs(make_key(k1, 1));
        SHUFFLE_TENSOR rhs(make_key(k2, 1));

        SHUFFLE_TENSOR expected;

        LET k12[] = {1, 2};
        LET k21[] = {2, 1};

        expected.add_scal_prod(make_key(k12, 2), 1.0);
        expected.add_scal_prod(make_key(k21, 2), 1.0);

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough, 5);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_and_add_deg_1_deg_1_op_max_depth

    TEST_FIXTURE(Fixture, test_multiply_deg_1_deg_1_op)
    {

        // test: {1(1)} * {1(2)} == {1(1,2) 1(2,1)}

        LET k1[] = {1};
        LET k2[] = {2};

        SHUFFLE_TENSOR lhs(make_key(k1, 1));
        SHUFFLE_TENSOR rhs(make_key(k2, 1));

        SHUFFLE_TENSOR expected;

        LET k12[] = {1, 2};
        LET k21[] = {2, 1};

        expected.add_scal_prod(make_key(k12, 2), 1.0);
        expected.add_scal_prod(make_key(k21, 2), 1.0);

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_deg_1_deg_1_op

    TEST_FIXTURE(Fixture, test_multiply_deg_1_deg_1_op_max_depth)
    {

        // test: {1(1)} * {1(2)} == {1(1,2) 1(2,1)}

        LET k1[] = {1};
        LET k2[] = {2};

        SHUFFLE_TENSOR lhs(make_key(k1, 1));
        SHUFFLE_TENSOR rhs(make_key(k2, 1));

        SHUFFLE_TENSOR expected;

        LET k12[] = {1, 2};
        LET k21[] = {2, 1};

        expected.add_scal_prod(make_key(k12, 2), 1.0);
        expected.add_scal_prod(make_key(k21, 2), 1.0);

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough, 5);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_deg_1_deg_1_op_max_depth

    TEST_FIXTURE(Fixture, test_multiply_inplace_deg_1_deg_1_op)
    {

        // test: {1(1)} * {1(2)} == {1(1,2) 1(2,1)}

        LET k1[] = {1};
        LET k2[] = {2};

        SHUFFLE_TENSOR lhs(make_key(k1, 1));
        SHUFFLE_TENSOR rhs(make_key(k2, 1));

        SHUFFLE_TENSOR expected;

        LET k12[] = {1, 2};
        LET k21[] = {2, 1};

        expected.add_scal_prod(make_key(k12, 2), 1.0);
        expected.add_scal_prod(make_key(k21, 2), 1.0);

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_inplace(my_shuffle_tensor_product, lhs, rhs, my_scalar_passthrough);

        CHECK_EQUAL(expected, lhs);

    }// TEST test_multiply_inplace_deg_1_deg_1_op

    TEST_FIXTURE(Fixture, test_multiply_inplace_deg_1_deg_1_op_max_depth)
    {

        // test: {1(1)} * {1(2)} == {1(1,2) 1(2,1)}

        LET k1[] = {1};
        LET k2[] = {2};

        SHUFFLE_TENSOR lhs(make_key(k1, 1));
        SHUFFLE_TENSOR rhs(make_key(k2, 1));

        SHUFFLE_TENSOR expected;

        LET k12[] = {1, 2};
        LET k21[] = {2, 1};

        expected.add_scal_prod(make_key(k12, 2), 1.0);
        expected.add_scal_prod(make_key(k21, 2), 1.0);

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_inplace(my_shuffle_tensor_product, lhs, rhs, my_scalar_passthrough, 5);

        CHECK_EQUAL(expected, lhs);

    }// TEST test_multiply_inplace_deg_1_deg_1_op_max_depth

    // test: {1(1)} * {1(2,3)} == {1(1,2,3) 1(2,1,3) 1(2,3,1)}

    TEST_FIXTURE(Fixture, test_star_deg_1_deg_2)
    {

        // test: {1(1)} * {1(2,3)} == {1(1,2,3) 1(2,1,3) 1(2,3,1)}

        LET k1[] = {1};
        LET k2[] = {2, 3};

        SHUFFLE_TENSOR lhs(make_key(k1, 1));// {1()1}
        SHUFFLE_TENSOR rhs(make_key(k2, 2));// {1(2,3)}

        SHUFFLE_TENSOR expected;

        LET k123[] = {1, 2, 3};
        LET k213[] = {2, 1, 3};
        LET k231[] = {2, 3, 1};

        expected.add_scal_prod(make_key(k123, 3), 1.0);
        expected.add_scal_prod(make_key(k213, 3), 1.0);
        expected.add_scal_prod(make_key(k231, 3), 1.0);

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;

        result = lhs * rhs;

        CHECK_EQUAL(expected, result);

    }// TEST test_star_unidim_deg_1_deg_2

    TEST_FIXTURE(Fixture, test_multiply_and_add_deg_1_deg_2_op)
    {

        // test: {1(1)} * {1(2,3)} == {1(1,2,3) 1(2,1,3) 1(2,3,1)}

        LET k1[] = {1};
        LET k2[] = {2, 3};

        SHUFFLE_TENSOR lhs(make_key(k1, 1));// {1()1}
        SHUFFLE_TENSOR rhs(make_key(k2, 2));// {1(2,3)}

        SHUFFLE_TENSOR expected;

        LET k123[] = {1, 2, 3};
        LET k213[] = {2, 1, 3};
        LET k231[] = {2, 3, 1};

        expected.add_scal_prod(make_key(k123, 3), 1.0);
        expected.add_scal_prod(make_key(k213, 3), 1.0);
        expected.add_scal_prod(make_key(k231, 3), 1.0);

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_and_add_deg_1_deg_2_op

    TEST_FIXTURE(Fixture, test_multiply_and_add_deg_1_deg_2_op_max_depth)
    {

        // test: {1(1)} * {1(2,3)} == {1(1,2,3) 1(2,1,3) 1(2,3,1)}

        LET k1[] = {1};
        LET k2[] = {2, 3};

        SHUFFLE_TENSOR lhs(make_key(k1, 1));// {1()1}
        SHUFFLE_TENSOR rhs(make_key(k2, 2));// {1(2,3)}

        SHUFFLE_TENSOR expected;

        LET k123[] = {1, 2, 3};
        LET k213[] = {2, 1, 3};
        LET k231[] = {2, 3, 1};

        expected.add_scal_prod(make_key(k123, 3), 1.0);
        expected.add_scal_prod(make_key(k213, 3), 1.0);
        expected.add_scal_prod(make_key(k231, 3), 1.0);

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough, 5);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_and_add_deg_1_deg_2_op_max_depth

    TEST_FIXTURE(Fixture, test_multiply_deg_1_deg_2_op)
    {

        // test: {1(1)} * {1(2,3)} == {1(1,2,3) 1(2,1,3) 1(2,3,1)}

        LET k1[] = {1};
        LET k2[] = {2, 3};

        SHUFFLE_TENSOR lhs(make_key(k1, 1));// {1()1}
        SHUFFLE_TENSOR rhs(make_key(k2, 2));// {1(2,3)}

        SHUFFLE_TENSOR expected;

        LET k123[] = {1, 2, 3};
        LET k213[] = {2, 1, 3};
        LET k231[] = {2, 3, 1};

        expected.add_scal_prod(make_key(k123, 3), 1.0);
        expected.add_scal_prod(make_key(k213, 3), 1.0);
        expected.add_scal_prod(make_key(k231, 3), 1.0);

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;


        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_deg_1_deg_2_op

    TEST_FIXTURE(Fixture, test_multiply_deg_1_deg_2_op_max_depth)
    {

        // test: {1(1)} * {1(2,3)} == {1(1,2,3) 1(2,1,3) 1(2,3,1)}

        LET k1[] = {1};
        LET k2[] = {2, 3};

        SHUFFLE_TENSOR lhs(make_key(k1, 1));// {1()1}
        SHUFFLE_TENSOR rhs(make_key(k2, 2));// {1(2,3)}

        SHUFFLE_TENSOR expected;

        LET k123[] = {1, 2, 3};
        LET k213[] = {2, 1, 3};
        LET k231[] = {2, 3, 1};

        expected.add_scal_prod(make_key(k123, 3), 1.0);
        expected.add_scal_prod(make_key(k213, 3), 1.0);
        expected.add_scal_prod(make_key(k231, 3), 1.0);

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough, 5);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_deg_1_deg_2_op_max_depth

    TEST_FIXTURE(Fixture, test_multiply_inplace_deg_1_deg_2_op)
    {

        // test: {1(1)} * {1(2,3)} == {1(1,2,3) 1(2,1,3) 1(2,3,1)}

        LET k1[] = {1};
        LET k2[] = {2, 3};

        SHUFFLE_TENSOR lhs(make_key(k1, 1));// {1()1}
        SHUFFLE_TENSOR rhs(make_key(k2, 2));// {1(2,3)}

        SHUFFLE_TENSOR expected;

        LET k123[] = {1, 2, 3};
        LET k213[] = {2, 1, 3};
        LET k231[] = {2, 3, 1};

        expected.add_scal_prod(make_key(k123, 3), 1.0);
        expected.add_scal_prod(make_key(k213, 3), 1.0);
        expected.add_scal_prod(make_key(k231, 3), 1.0);

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_inplace(my_shuffle_tensor_product, lhs, rhs, my_scalar_passthrough);

        CHECK_EQUAL(expected, lhs);

    }// TEST test_multiply_inplace_deg_1_deg_2_op

    TEST_FIXTURE(Fixture, test_multiply_inplace_deg_1_deg_2_op_max_depth)
    {

        // test: {1(1)} * {1(2,3)} == {1(1,2,3) 1(2,1,3) 1(2,3,1)}

        LET k1[] = {1};
        LET k2[] = {2, 3};

        SHUFFLE_TENSOR lhs(make_key(k1, 1));// {1()1}
        SHUFFLE_TENSOR rhs(make_key(k2, 2));// {1(2,3)}

        SHUFFLE_TENSOR expected;

        LET k123[] = {1, 2, 3};
        LET k213[] = {2, 1, 3};
        LET k231[] = {2, 3, 1};

        expected.add_scal_prod(make_key(k123, 3), 1.0);
        expected.add_scal_prod(make_key(k213, 3), 1.0);
        expected.add_scal_prod(make_key(k231, 3), 1.0);

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_inplace(my_shuffle_tensor_product, lhs, rhs, my_scalar_passthrough, 5);

        CHECK_EQUAL(expected, lhs);

    }// TEST test_multiply_inplace_deg_1_deg_2_op_max_depth

    // test: {1(1,2)} * {1(3)} == {1(1,2,3) 1(1,3,2) 1(3,1,2)}

    TEST_FIXTURE(Fixture, test_star_deg_2_deg_1)
    {

        // test: {1(1,2)} * {1(3)} == {1(1,2,3) 1(1,3,2) 1(3,1,2)}

        LET k1[] = {1, 2};
        LET k2[] = {3};

        SHUFFLE_TENSOR lhs(make_key(k1, 2));// {1(1,2)}
        SHUFFLE_TENSOR rhs(make_key(k2, 1));// {1(3)}

        SHUFFLE_TENSOR expected;

        LET k123[] = {1, 2, 3};
        LET k132[] = {1, 3, 2};
        LET k312[] = {3, 1, 2};

        expected.add_scal_prod(make_key(k123, 3), 1.0);// 1(1,2,3)
        expected.add_scal_prod(make_key(k132, 3), 1.0);// 1(1,3,2)
        expected.add_scal_prod(make_key(k312, 3), 1.0);// 1(3,1,2)      // {1(1,2,3) 1(1,3,2) 1(3,1,2)}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;

        result = lhs * rhs;

        CHECK_EQUAL(expected, result);

    }// TEST test_product_unidim_deg_2_deg_1

    TEST_FIXTURE(Fixture, test_multiply_and_add_deg_2_deg_1_op)
    {

        // test: {1(1,2)} * {1(3)} == {1(1,2,3) 1(1,3,2) 1(3,1,2)}

        LET k1[] = {1, 2};
        LET k2[] = {3};

        SHUFFLE_TENSOR lhs(make_key(k1, 2));// {1(1,2)}
        SHUFFLE_TENSOR rhs(make_key(k2, 1));// {1(3)}

        SHUFFLE_TENSOR expected;

        LET k123[] = {1, 2, 3};
        LET k132[] = {1, 3, 2};
        LET k312[] = {3, 1, 2};

        expected.add_scal_prod(make_key(k123, 3), 1.0);// 1(1,2,3)
        expected.add_scal_prod(make_key(k132, 3), 1.0);// 1(1,3,2)
        expected.add_scal_prod(make_key(k312, 3), 1.0);// 1(3,1,2)      // {1(1,2,3) 1(1,3,2) 1(3,1,2)}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_and_add_deg_2_deg_1_op

    TEST_FIXTURE(Fixture, test_multiply_and_add_deg_2_deg_1_op_max_depth)
    {

        // test: {1(1,2)} * {1(3)} == {1(1,2,3) 1(1,3,2) 1(3,1,2)}

        LET k1[] = {1, 2};
        LET k2[] = {3};

        SHUFFLE_TENSOR lhs(make_key(k1, 2));// {1(1,2)}
        SHUFFLE_TENSOR rhs(make_key(k2, 1));// {1(3)}

        SHUFFLE_TENSOR expected;

        LET k123[] = {1, 2, 3};
        LET k132[] = {1, 3, 2};
        LET k312[] = {3, 1, 2};

        expected.add_scal_prod(make_key(k123, 3), 1.0);// 1(1,2,3)
        expected.add_scal_prod(make_key(k132, 3), 1.0);// 1(1,3,2)
        expected.add_scal_prod(make_key(k312, 3), 1.0);// 1(3,1,2)      // {1(1,2,3) 1(1,3,2) 1(3,1,2)}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough, 5);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_and_add_deg_2_deg_1_op_max_depth

    TEST_FIXTURE(Fixture, test_multiply_deg_2_deg_1_op)
    {

        // test: {1(1,2)} * {1(3)} == {1(1,2,3) 1(1,3,2) 1(3,1,2)}

        LET k1[] = {1, 2};
        LET k2[] = {3};

        SHUFFLE_TENSOR lhs(make_key(k1, 2));// {1(1,2)}
        SHUFFLE_TENSOR rhs(make_key(k2, 1));// {1(3)}

        SHUFFLE_TENSOR expected;

        LET k123[] = {1, 2, 3};
        LET k132[] = {1, 3, 2};
        LET k312[] = {3, 1, 2};

        expected.add_scal_prod(make_key(k123, 3), 1.0);// 1(1,2,3)
        expected.add_scal_prod(make_key(k132, 3), 1.0);// 1(1,3,2)
        expected.add_scal_prod(make_key(k312, 3), 1.0);// 1(3,1,2)      // {1(1,2,3) 1(1,3,2) 1(3,1,2)}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_deg_2_deg_1_op

    TEST_FIXTURE(Fixture, test_multiply_deg_2_deg_1_op_max_depth)
    {

        // test: {1(1,2)} * {1(3)} == {1(1,2,3) 1(1,3,2) 1(3,1,2)}

        LET k1[] = {1, 2};
        LET k2[] = {3};

        SHUFFLE_TENSOR lhs(make_key(k1, 2));// {1(1,2)}
        SHUFFLE_TENSOR rhs(make_key(k2, 1));// {1(3)}

        SHUFFLE_TENSOR expected;

        LET k123[] = {1, 2, 3};
        LET k132[] = {1, 3, 2};
        LET k312[] = {3, 1, 2};

        expected.add_scal_prod(make_key(k123, 3), 1.0);// 1(1,2,3)
        expected.add_scal_prod(make_key(k132, 3), 1.0);// 1(1,3,2)
        expected.add_scal_prod(make_key(k312, 3), 1.0);// 1(3,1,2)      // {1(1,2,3) 1(1,3,2) 1(3,1,2)}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR result;
        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_and_add(my_shuffle_tensor_product, result, lhs, rhs, my_scalar_passthrough, 5);

        CHECK_EQUAL(expected, result);

    }// TEST test_multiply_deg_2_deg_1_op_max_depth

    TEST_FIXTURE(Fixture, test_multiply_inplace_deg_2_deg_1_op)
    {

        // test: {1(1,2)} * {1(3)} == {1(1,2,3) 1(1,3,2) 1(3,1,2)}

        LET k1[] = {1, 2};
        LET k2[] = {3};

        SHUFFLE_TENSOR lhs(make_key(k1, 2));// {1(1,2)}
        SHUFFLE_TENSOR rhs(make_key(k2, 1));// {1(3)}

        SHUFFLE_TENSOR expected;

        LET k123[] = {1, 2, 3};
        LET k132[] = {1, 3, 2};
        LET k312[] = {3, 1, 2};

        expected.add_scal_prod(make_key(k123, 3), 1.0);// 1(1,2,3)
        expected.add_scal_prod(make_key(k132, 3), 1.0);// 1(1,3,2)
        expected.add_scal_prod(make_key(k312, 3), 1.0);// 1(3,1,2)      // {1(1,2,3) 1(1,3,2) 1(3,1,2)}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_inplace(my_shuffle_tensor_product, lhs, rhs, my_scalar_passthrough);

        CHECK_EQUAL(expected, lhs);

    }// TEST test_multiply_inplace_deg_2_deg_1_op

    TEST_FIXTURE(Fixture, test_multiply_inplace_deg_2_deg_1_op_max_depth)
    {

        // test: {1(1,2)} * {1(3)} == {1(1,2,3) 1(1,3,2) 1(3,1,2)}

        LET k1[] = {1, 2};
        LET k2[] = {3};

        SHUFFLE_TENSOR lhs(make_key(k1, 2));// {1(1,2)}
        SHUFFLE_TENSOR rhs(make_key(k2, 1));// {1(3)}

        SHUFFLE_TENSOR expected;

        LET k123[] = {1, 2, 3};
        LET k132[] = {1, 3, 2};
        LET k312[] = {3, 1, 2};

        expected.add_scal_prod(make_key(k123, 3), 1.0);// 1(1,2,3)
        expected.add_scal_prod(make_key(k132, 3), 1.0);// 1(1,3,2)
        expected.add_scal_prod(make_key(k312, 3), 1.0);// 1(3,1,2)      // {1(1,2,3) 1(1,3,2) 1(3,1,2)}

        // std::cout << "lhs = " << lhs << ", rhs = " << rhs << ", expected = " << expected << std::endl;

        SHUFFLE_TENSOR_MULTIPLICATION my_shuffle_tensor_product;
        alg::mult::scalar_passthrough my_scalar_passthrough;
        mtraits::multiply_inplace(my_shuffle_tensor_product, lhs, rhs, my_scalar_passthrough, 5);

        CHECK_EQUAL(expected, lhs);

    }// TEST test_multiply_inplace_deg_2_deg_1_op_max_depth

    // test: < Sh(t1), S(x) > x < Sh(t2), S(x) > = <Sh(t1) x Sh(t2), S(x) >

    TEST_FIXTURE(Fixture, test_shuffle_product_accuracy)
    {

        LET k1[] = {1};
        LET k2[] = {2};
        LET k3[] = {3};
        LET k4[] = {4};
        LET k5[] = {5};

        TENSOR t1(make_key(k1, 1));
        TENSOR t2(make_key(k2, 1));
        TENSOR t3(make_key(k3, 1));
        TENSOR t4(make_key(k4, 1));
        TENSOR t5(make_key(k5, 1));

        // std::cout << "t1 = " << t1 << ", t2 = " << t2 << ", t3 = " << t3 << ", t4 = " << t4 << ", t5 = " << t5 << std::endl;

        TENSOR sig = exp(t1) * exp(t2) * exp(t3) * exp(t4) * exp(t5);

        //std::cout << "signature = " << sig << std::endl;

        LET k11[] = {1, 1};
        LET k111[] = {1, 1, 1};

        SHUFFLE_TENSOR st1;

        st1.add_scal_prod(make_key(k11, 2), 1);
        st1.add_scal_prod(make_key(k111, 3), 1);

        LET k22[] = {2, 2};
        //LET k222[] = {2, 2, 2};

        SHUFFLE_TENSOR st2;

        st2.add_scal_prod(make_key(k22, 2), 1);
        //st2.add_scal_prod(make_key(k222, 3), 1);

        //std::cout << "st1 = " << st1 << "\nst2 = " << st2 << "\nst1 *st2 = " << st1 * st2 << std::endl;

        pairing<COEFF, 5, 5> my_pairing;

        COEFF::S lhs;
        COEFF::S rhs;

        /*
         * Something interesting happens here. The shuffle tensor truncates the terms of degree 6, which leads to
         * a non-zero error term in the inner product between <sh1, sig><sh2, sig> and <sh1*sh2, sig>. I've trimmed sh2
         * so that it only has degree 2 terms so no truncation occurs, and we should get a precise answer.
         */

        lhs = my_pairing(st1, sig) * my_pairing(st2, sig);
        rhs = my_pairing(st1 * st2, sig);

        CHECK_EQUAL(lhs, rhs);

    }// test_shuffle_product_accuracy

}// SUITE shuffle_tensor
