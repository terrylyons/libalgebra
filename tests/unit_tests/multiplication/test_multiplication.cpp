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

    static constexpr DEG WIDTH1 = WIDTH; // note added 1 to name, error: declaration of 'constexpr const DEG WIDTH' shadows template parameter
    static constexpr DEG DEPTH1 = DEPTH; // note added 1 to name, error: declaration of 'constexpr const DEG DEPTH' shadows template parameter
    static constexpr DEG poly_width = hall_basis<WIDTH, DEPTH>::start_of_degree(DEPTH + 1);

    using poly_t = alg::poly<scalar_field>;
    using poly_coeffs = coefficients::coefficient_ring<poly_t, typename scalar_field::Q>;
    static_assert(std::is_same<poly_coeffs::S, poly_t>::value, "the trait class of a coefficient ring must report the type of the coefficients");

    template<DEG DEPTH1> // note added 1 to name, error: declaration of 'constexpr const DEG DEPTH' shadows template parameter
    using LIE_ = alg::lie<poly_coeffs, WIDTH, DEPTH, vectors::dense_vector>;
    using LIE = LIE_<DEPTH>;

    template<DEG DEPTH1> // note added 1 to name, error: declaration of 'constexpr const DEG DEPTH' shadows template parameter
    using lie_basis_t_ = lie_basis<WIDTH, DEPTH>;

    using lie_basis_t = lie_basis_t_<WIDTH>;
    lie_basis_t lbasis;

    template<DEG DEPTH1> // note added 1 to name, error: declaration of 'constexpr const DEG DEPTH' shadows template parameter
    using TENSOR_ = alg::free_tensor<poly_coeffs, WIDTH, DEPTH, vectors::dense_vector>;
    using TENSOR = TENSOR_<DEPTH>;

    template<DEG DEPTH1> // note added 1 to name, error: declaration of 'constexpr const DEG DEPTH' shadows template parameter
    using tensor_basis_t_ = alg::tensor_basis<WIDTH, DEPTH>;

    using tensor_basis_t = tensor_basis_t_< DEPTH>;
    tensor_basis_t tbasis;

    //using SHUFFLE_TENSOR = alg::shuffle_tensor<scalar_field, WIDTH, DEPTH>;
    template<DEG DEPTH1> // note added 1 to name, error: declaration of 'constexpr const DEG DEPTH' shadows template parameter
    using SHUFFLE_TENSOR_ = alg::shuffle_tensor<poly_coeffs, WIDTH, DEPTH>;
    using SHUFFLE_TENSOR = SHUFFLE_TENSOR_<DEPTH>;

    template<DEG DEPTH1> // note added 1 to name, error: declaration of 'constexpr const DEG DEPTH' shadows template parameter
    using shuffle_tensor_basis_t_ = alg::tensor_basis<WIDTH, DEPTH>;
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
//        if constexpr (std::is_integral<decltype(basis.begin())>::value) {
//            // invalid code if the premise is false - constexpr essential to avoid compilation
//            count = basis.begin();
//        }
//        else {
//            count = 0;
//        } // TODO: Fix

        std::map<int, std::pair<typename VECTOR::KEY, std::string>> legend;

        VECTOR result;
        //for range type for loop approach use ": basis.iterate_keys()"
        for (auto key = basis.begin(), end = basis.end(); key != end; key = basis.nextkey(key)) {
            result[key] = poly_t(count + offset, 1);

            // record the mapping from keys to monomials
            auto basis_key_pair = std::pair<typename VECTOR::BASIS*, typename VECTOR::KEY>(&VECTOR::basis, key);
            std::stringstream buffer;
            buffer << basis_key_pair;
//            legend[count + offset] = std::pair(key, buffer.str());  // TODO: Fix
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

    TEST_FIXTURE(dense_fixture, cmake_check)
    {
        CHECK_EQUAL(0, 1);
    }

    // the receiving environment
    constexpr DEG WIDTHOUT = 3;
    constexpr DEG DEPTHOUT = 2;

// the depth of the shuffles producing the channels
    constexpr DEG INOUTDEPTH = 2;

// the incoming stream with enough accuracy to determine the
// required integrals of the projections
    constexpr DEG WIDTHIN = 2;
    constexpr DEG DEPTHIN = (DEPTHOUT * INOUTDEPTH);

// the input and output environments
    using IN = Environment<WIDTHIN, DEPTHIN>;
    using OUT = Environment<WIDTHOUT, DEPTHOUT>;

    // the input and output environments
    using IN = Environment<WIDTHIN, DEPTHIN>;
    using OUT = Environment<WIDTHOUT, DEPTHOUT>;

    // steady state requires inputs limited to the same depth as the output
    using SHORT_LIE = IN::LIE_<DEPTHOUT>;

// limiting the degree of nonlinearity in the functions on paths
    using SHORT_SHUFFLE = IN::SHUFFLE_TENSOR_<INOUTDEPTH>;

//    using SHORT_LIE = environment_fixture::LIE_<4>;

    TEST_FIXTURE(IN, environment_check)
    {
        IN in;
        const auto& sbasis = in.sbasis;

        std::cout << "Creating the two generic input log signatures \"before\" and \"during\" truncated to level " << DEPTHIN << "\n\n";
        IN::LIE logsig_before;
        add_equals_short(logsig_before, in.generic_vector<SHORT_LIE>(1000));
        SHOW(logsig_before);
        IN::TENSOR tensor_logsig_before = in.maps_.l2t(logsig_before);
        IN::TENSOR sig_before = exp(tensor_logsig_before);
        SHOW(sig_before);

//        in.generic_vector<SHORT_LIE>(1000);
//
//        SHOW(in.generic_vector<SHORT_LIE>(1000));



        CHECK_EQUAL(0, 1);
    }

}// SUITE Multiplication