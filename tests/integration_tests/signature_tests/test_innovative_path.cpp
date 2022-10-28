//
// Created by sam on 11/03/2021.
//

#include <algorithm>
#include <map>
#include <vector>

#include <UnitTest++.h>

#include "libalgebra/alg_types.h"
#include "libalgebra/libalgebra.h"

#include "tests/common/helpers.h"
#include "tests/common/rng.h"
#include "tests/common/time_and_details.h"

#include "generic_coefficient.h"
#include "generic_lie_increment.h"
#include "generic_path.h"

template<unsigned Width>
generic_path<Width> make_innovative_path(unsigned length)
{
    std::vector<generic_lie_increment<Width, int32_t>> tmp;
    tmp.reserve(length);
    mt19937 rng;
    UNIFORM_INT_DIST<int32_t> dist(0, Width - 1);

    UNIFORM_INT_DIST<int32_t> distn(-5, 5);
    UNIFORM_INT_DIST<int32_t> distd(1, 60);

    std::vector<size_t> seen;
    seen.reserve(length);

    size_t key = dist(rng);
    std::vector<generic_coefficient<int32_t>> incr_tmp;
    incr_tmp.reserve(Width);
    for (size_t i = 0; i < Width; ++i) {
        if (i == key) {
            incr_tmp.push_back(generic_coefficient<int32_t>(distn(rng), distd(rng)));
        }
        else {
            incr_tmp.push_back(generic_coefficient<int32_t>(0, 1));
        }
    }
    tmp.push_back(generic_lie_increment<Width, int32_t>(incr_tmp));
    seen.push_back(key);

    for (size_t i = 1; i < length; ++i) {
        while (std::find(seen.begin(), seen.end(), key) != seen.end()) {
            key = dist(rng);
        }
        incr_tmp.clear();
        incr_tmp.reserve(Width);
        for (size_t i = 0; i < Width; ++i) {
            if (i == key) {
                incr_tmp.push_back(generic_coefficient<int32_t>(distn(rng), distd(rng)));
            }
            else {
                incr_tmp.push_back(generic_coefficient<int32_t>(0, 1));
            }
        }
        tmp.push_back(generic_lie_increment<Width, int32_t>(incr_tmp));
        ;
        seen.push_back(key);
    }
    assert(tmp.size() == length);

    return generic_path<Width>(tmp);
}

namespace {

#if __cplusplus >= 201103UL
constexpr
#endif
        inline size_t
        binomial_coeff(size_t n, size_t k)
{
    return (k > n) ? 0 : (k == 0 || k == n) ? 1
                                            : binomial_coeff(n - 1, k - 1) + binomial_coeff(n - 1, k);
}

}// namespace

template<unsigned Width, unsigned Depth>
struct GenericFixture {
    static const unsigned width = Width;
    static const unsigned depth = Depth;
    generic_path<width> path;

    double expected_double_error;
    float expected_float_error;

    typedef typename alg_types<2, 2, Rational>::SCA Rat;

    typedef alg::coefficients::coefficient_field<Rat> rational_field;
    typedef alg::coefficients::double_field double_field;
    typedef alg::coefficients::float_field float_field;

    struct rational_sparse_framework {
        typedef typename rational_field::S S;
        typedef typename rational_field::Q Q;
        typedef alg::free_tensor_basis<width, depth> TBASIS;
        typedef alg::lie_basis<width, depth> LBASIS;

        typedef alg::vectors::sparse_vector<
                TBASIS,
                rational_field,
                std::map<typename TBASIS::KEY, S>>
                SPTENS;

        typedef alg::vectors::sparse_vector<
                LBASIS,
                rational_field,
                std::map<typename LBASIS::KEY, S>>
                SPLIE;

        typedef alg::free_tensor<rational_field, width, depth, alg::vectors::sparse_vector> TENSOR;
        typedef alg::lie<rational_field, width, depth, alg::vectors::sparse_vector> LIE;
        typedef alg::maps<rational_field, width, depth, TENSOR, LIE> MAPS;
    };

    struct rational_dense_framework {
        typedef typename rational_field::S S;
        typedef typename rational_field::Q Q;
        typedef alg::free_tensor_basis<width, depth> TBASIS;
        typedef alg::lie_basis<width, depth> LBASIS;

        typedef alg::vectors::dense_vector<
                TBASIS,
                rational_field>
                DTENS;

        typedef alg::vectors::sparse_vector<
                LBASIS,
                rational_field>
                DLIE;

        typedef alg::free_tensor<rational_field, width, depth, alg::vectors::dense_vector> TENSOR;
        typedef alg::lie<rational_field, width, depth, alg::vectors::dense_vector> LIE;
        typedef alg::maps<rational_field, width, depth, TENSOR, LIE> MAPS;
    };

    struct rational_hybrid_framework {
        typedef typename rational_field::S S;
        typedef typename rational_field::Q Q;
        typedef alg::free_tensor_basis<width, depth> TBASIS;
        typedef alg::lie_basis<width, depth> LBASIS;

        typedef alg::vectors::hybrid_vector<
                TBASIS,
                rational_field,
                alg::vectors::policy::basic_resize_policy,
                std::map<typename TBASIS::KEY, S>>
                HTENS;

        typedef alg::vectors::hybrid_vector<
                LBASIS,
                rational_field,
                alg::vectors::policy::basic_resize_policy,
                std::map<typename LBASIS::KEY, S>>
                HLIE;

        typedef alg::free_tensor<rational_field, width, depth, alg::vectors::hybrid_vector> TENSOR;
        typedef alg::lie<rational_field, width, depth, alg::vectors::hybrid_vector> LIE;
        typedef alg::maps<rational_field, width, depth, TENSOR, LIE> MAPS;
    };

    struct double_sparse_framework {
        typedef typename double_field::S S;
        typedef typename double_field::Q Q;
        typedef alg::free_tensor_basis<width, depth> TBASIS;
        typedef alg::lie_basis<width, depth> LBASIS;

        typedef alg::vectors::sparse_vector<
                TBASIS,
                double_field,
                std::map<typename TBASIS::KEY, S>>
                SPTENS;

        typedef alg::vectors::sparse_vector<
                LBASIS,
                double_field,
                std::map<typename LBASIS::KEY, S>>
                SPLIE;

        typedef alg::free_tensor<double_field, width, depth, alg::vectors::sparse_vector> TENSOR;
        typedef alg::lie<double_field, width, depth, alg::vectors::sparse_vector> LIE;
        typedef alg::maps<double_field, width, depth, TENSOR, LIE> MAPS;
    };

    struct double_dense_framework {
        typedef typename double_field::S S;
        typedef typename double_field::Q Q;
        typedef alg::free_tensor_basis<width, depth> TBASIS;
        typedef alg::lie_basis<width, depth> LBASIS;

        typedef alg::vectors::dense_vector<
                TBASIS,
                double_field>
                DTENS;

        typedef alg::vectors::sparse_vector<
                LBASIS,
                double_field>
                DLIE;

        typedef alg::free_tensor<double_field, width, depth, alg::vectors::dense_vector> TENSOR;
        typedef alg::lie<double_field, width, depth, alg::vectors::dense_vector> LIE;
        typedef alg::maps<double_field, width, depth, TENSOR, LIE> MAPS;
    };

    struct double_hybrid_framework {
        typedef typename double_field::S S;
        typedef typename double_field::Q Q;
        typedef alg::free_tensor_basis<width, depth> TBASIS;
        typedef alg::lie_basis<width, depth> LBASIS;

        typedef alg::vectors::hybrid_vector<
                TBASIS,
                double_field,
                alg::vectors::policy::basic_resize_policy,
                std::map<typename TBASIS::KEY, S>>
                HTENS;

        typedef alg::vectors::hybrid_vector<
                LBASIS,
                double_field,
                alg::vectors::policy::basic_resize_policy,
                std::map<typename LBASIS::KEY, S>>
                HLIE;

        typedef alg::free_tensor<double_field, width, depth, alg::vectors::hybrid_vector> TENSOR;
        typedef alg::lie<double_field, width, depth, alg::vectors::hybrid_vector> LIE;
        typedef alg::maps<double_field, width, depth, TENSOR, LIE> MAPS;
    };

    struct float_sparse_framework {
        typedef typename float_field::S S;
        typedef typename float_field::Q Q;
        typedef alg::free_tensor_basis<width, depth> TBASIS;
        typedef alg::lie_basis<width, depth> LBASIS;

        typedef alg::vectors::sparse_vector<
                TBASIS,
                float_field,
                std::map<typename TBASIS::KEY, S>>
                SPTENS;

        typedef alg::vectors::sparse_vector<
                LBASIS,
                float_field,
                std::map<typename LBASIS::KEY, S>>
                SPLIE;

        typedef alg::free_tensor<float_field, width, depth, alg::vectors::sparse_vector> TENSOR;
        typedef alg::lie<float_field, width, depth, alg::vectors::sparse_vector> LIE;
        typedef alg::maps<float_field, width, depth, TENSOR, LIE> MAPS;
    };

    struct float_dense_framework {
        typedef typename float_field::S S;
        typedef typename float_field::Q Q;
        typedef alg::free_tensor_basis<width, depth> TBASIS;
        typedef alg::lie_basis<width, depth> LBASIS;

        typedef alg::vectors::dense_vector<
                TBASIS,
                float_field>
                DTENS;

        typedef alg::vectors::sparse_vector<
                LBASIS,
                float_field>
                DLIE;

        typedef alg::free_tensor<float_field, width, depth, alg::vectors::dense_vector> TENSOR;
        typedef alg::lie<float_field, width, depth, alg::vectors::dense_vector> LIE;
        typedef alg::maps<float_field, width, depth, TENSOR, LIE> MAPS;
    };

    struct float_hybrid_framework {
        typedef typename float_field::S S;
        typedef typename float_field::Q Q;
        typedef alg::free_tensor_basis<width, depth> TBASIS;
        typedef alg::lie_basis<width, depth> LBASIS;

        typedef alg::vectors::hybrid_vector<
                TBASIS,
                float_field,
                alg::vectors::policy::basic_resize_policy,
                std::map<typename TBASIS::KEY, S>>
                HTENS;

        typedef alg::vectors::hybrid_vector<
                LBASIS,
                float_field,
                alg::vectors::policy::basic_resize_policy,
                std::map<typename LBASIS::KEY, S>>
                HLIE;

        typedef alg::free_tensor<float_field, width, depth, alg::vectors::hybrid_vector> TENSOR;
        typedef alg::lie<float_field, width, depth, alg::vectors::hybrid_vector> LIE;
        typedef alg::maps<float_field, width, depth, TENSOR, LIE> MAPS;
    };

    GenericFixture() : path(make_innovative_path<width>(width)),
                       expected_double_error(1.0e-14),
                       expected_float_error(1.0e-6f)
    {}

    size_t sig_support(unsigned increments = width)
    {
        return binomial_coeff(increments + depth, depth);
    }
};

SUITE(innovative_path_5_tests)
{

    typedef GenericFixture<5, 5> Fixture;

#include "accuracy_suite.ins"
#include "double_path_suite.ins"
#include "float_path_suite.ins"
#include "rational_path_suite.ins"
}

SUITE(innovative_path_10_tests)
{

    typedef GenericFixture<10, 3> Fixture;

#include "accuracy_suite.ins"
#include "double_path_suite.ins"
#include "float_path_suite.ins"
#include "rational_path_suite.ins"
}

SUITE(innovative_path_15_tests)
{

    typedef GenericFixture<15, 2> Fixture;

#include "accuracy_suite.ins"
#include "double_path_suite.ins"
#include "float_path_suite.ins"
#include "rational_path_suite.ins"
}
