//
// Created by sam on 16/03/2021.
//

#include <algorithm>
#include <cmath>
#include <exception>
#include <map>
#include <sstream>
#include <utility>
#include <vector>

#pragma warning(suppress : 4616)
#include <UnitTest++.h>
#include <boost/integer/common_factor_rt.hpp>// deprecated
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/special_functions/modf.hpp>

#include "libalgebra/alg_types.h"
#include "libalgebra/libalgebra.h"

#include "tests/common/helpers.h"
#include "tests/common/memfile.h"
#include "tests/common/rng.h"
#include "tests/common/time_and_details.h"

#include "generic_coefficient.h"
#include "generic_lie_increment.h"
#include "generic_path.h"

namespace {

template<typename Integer>
generic_coefficient<Integer> convert_to_generic(double arg, Integer prec)
{
    double integ = NAN;
    double fract = boost::math::modf(arg, &integ);

    assert(!boost::math::isnan(integ));

    Integer integral_part = static_cast<Integer>(integ);
    Integer scaled_fract = static_cast<Integer>(round(fract * static_cast<double>(prec)));

    Integer gcd1 = boost::integer::gcd(scaled_fract, prec);
    assert(gcd1 != 0);

    Integer fract_part_num = scaled_fract / gcd1;
    Integer fract_part_den = prec / gcd1;

    Integer tmp = fract_part_den * integral_part + fract_part_num;
    Integer gcd2 = boost::integer::gcd(tmp, fract_part_den);
    assert(gcd2 != 0);
    return generic_coefficient<Integer>(tmp / gcd2, fract_part_den / gcd2);
}

}// namespace

template<unsigned Width, typename Integer, typename Rng>
inline generic_lie_increment<Width, Integer>
make_brownian_increment(double stdev, Rng rng)
{
    std::vector<generic_coefficient<Integer>> tmp;
    tmp.reserve(Width);

    NORMAL_DIST<double> normal(0.0, stdev);

    for (unsigned i = 0; i < Width; ++i) {
        tmp.push_back(convert_to_generic<Integer>(normal(rng), 1000));
    }

    return generic_lie_increment<Width, Integer>(tmp);
}

template<unsigned Width>
generic_path<Width> make_brownian_path(size_t length = Width)
{
    typedef generic_lie_increment<Width, int32_t> value_type;
    assert(length > 0);
    static std::vector<value_type> tmp;
    if (tmp.empty()) {
        tmp.reserve(length);

        mt19937 rng(12345);
        double stdev = 1.0 / sqrt(static_cast<double>(length));

        for (std::size_t i = 0; i < length; ++i) {
            value_type incr(make_brownian_increment<Width, int32_t>(stdev, rng));
            tmp.push_back(incr);
        }
    }

    assert(tmp.size() == length);

    return generic_path<Width>(tmp);
}

template<unsigned Width, unsigned Depth, unsigned Length>
struct GenericFixture {
    static const unsigned width = Width;
    static const unsigned depth = Depth;
    static const unsigned length = Length;

    generic_path<Width> path;

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

    GenericFixture()
    {
    }

    size_t sig_support(unsigned increments = length)
    {
        if (increments == 0) {
            return 1;
        }
        else {
            typename rational_dense_framework::TBASIS basis;
            return basis.start_of_degree(depth + 1);
        }
    }
};

SUITE(brownian_path_5_5_10_tests)
{

    static const float expected_float_error = 4e-5f;
    static const double expected_double_error = 2e-12;

    struct Fixture : GenericFixture<5, 5, 10> {
        typedef GenericFixture<5, 5, 10> base;

        using base::rational_dense_framework;
        using base::rational_field;
        using base::rational_hybrid_framework;
        using base::rational_sparse_framework;

        using base::double_dense_framework;
        using base::double_field;
        using base::double_hybrid_framework;
        using base::double_sparse_framework;

        using base::float_dense_framework;
        using base::float_field;
        using base::float_hybrid_framework;
        using base::float_sparse_framework;

        Fixture()
        {
            int32_t vals[][5][2] = {
                    {{81, 200}, {39, 125}, {-43, 1000}, {-71, 500}, {76, 125}},
                    {{3, 250}, {61, 125}, {-71, 250}, {-29, 100}, {13, 500}},
                    {{-313, 500}, {187, 500}, {73, 250}, {-43, 250}, {127, 1000}},
                    {{-379, 500}, {17, 200}, {-43, 200}, {39, 200}, {17, 125}},
                    {{459, 1000}, {-207, 500}, {-117, 200}, {-123, 250}, {-149, 1000}},
                    {{-297, 500}, {-27, 50}, {3, 50}, {-69, 200}, {-11, 40}},
                    {{-291, 1000}, {63, 1000}, {-3, 50}, {-203, 1000}, {71, 200}},
                    {{21, 500}, {99, 500}, {11, 1000}, {67, 1000}, {31, 500}},
                    {{39, 1000}, {2, 125}, {-81, 250}, {-7, 250}, {-227, 500}},
                    {{19, 40}, {-53, 250}, {-1, 100}, {179, 1000}, {41, 200}}};

            path = generic_path<5>(reinterpret_cast<int32_t*>(vals), 10);
        }
    };

#include "accuracy_suite.ins"
#include "double_path_suite.ins"
#include "float_path_suite.ins"
#include "rational_path_suite.ins"
}

SUITE(brownian_path_5_5_50_tests)
{

    static const float expected_float_error = 7e-3f;
    static const double expected_double_error = 2e-11;

    struct Fixture : GenericFixture<5, 5, 50> {
        typedef GenericFixture<5, 5, 50> base;

        using base::rational_dense_framework;
        using base::rational_field;
        using base::rational_hybrid_framework;
        using base::rational_sparse_framework;

        using base::double_dense_framework;
        using base::double_field;
        using base::double_hybrid_framework;
        using base::double_sparse_framework;

        using base::float_dense_framework;
        using base::float_field;
        using base::float_hybrid_framework;
        using base::float_sparse_framework;

        Fixture()
        {
            int32_t vals[][5][2] = {
                    {{-43, 500}, {67, 500}, {-303, 1000}, {43, 500}, {-33, 1000}},
                    {{81, 500}, {-3, 250}, {2, 125}, {-133, 1000}, {149, 500}},
                    {{-93, 1000}, {7, 200}, {237, 1000}, {21, 200}, {137, 1000}},
                    {{-117, 1000}, {109, 1000}, {11, 1000}, {47, 1000}, {0, 1}},
                    {{-113, 1000}, {-59, 1000}, {-51, 1000}, {3, 50}, {-153, 1000}},
                    {{59, 1000}, {-9, 50}, {-29, 500}, {13, 250}, {-111, 1000}},
                    {{29, 250}, {-19, 500}, {-221, 1000}, {-89, 500}, {-51, 500}},
                    {{-53, 250}, {77, 1000}, {27, 1000}, {-7, 200}, {-83, 500}},
                    {{-47, 500}, {-41, 250}, {-177, 1000}, {67, 1000}, {81, 1000}},
                    {{-173, 1000}, {-1, 5}, {-1, 50}, {-43, 250}, {-21, 200}},
                    {{131, 1000}, {39, 500}, {123, 500}, {27, 200}, {-23, 100}},
                    {{69, 500}, {-7, 200}, {11, 125}, {-261, 1000}, {-3, 1000}},
                    {{87, 1000}, {-4, 125}, {211, 1000}, {-11, 125}, {79, 1000}},
                    {{-41, 500}, {-1, 20}, {23, 500}, {9, 100}, {8, 125}},
                    {{-87, 500}, {71, 500}, {1, 10}, {-13, 200}, {59, 1000}},
                    {{-2, 25}, {-113, 1000}, {-37, 250}, {-31, 500}, {77, 200}},
                    {{29, 500}, {91, 1000}, {-101, 1000}, {9, 125}, {-89, 1000}},
                    {{23, 250}, {-101, 1000}, {51, 500}, {33, 200}, {277, 1000}},
                    {{-17, 250}, {-81, 1000}, {53, 1000}, {51, 500}, {87, 1000}},
                    {{-27, 500}, {-13, 100}, {143, 500}, {27, 200}, {93, 500}},
                    {{-41, 1000}, {-77, 250}, {61, 250}, {3, 100}, {67, 500}},
                    {{4, 125}, {91, 1000}, {47, 1000}, {-6, 125}, {109, 500}},
                    {{87, 1000}, {-13, 1000}, {13, 125}, {-47, 250}, {-57, 500}},
                    {{-7, 1000}, {-71, 500}, {11, 250}, {-73, 1000}, {11, 500}},
                    {{71, 1000}, {-11, 1000}, {6, 125}, {-3, 125}, {21, 200}},
                    {{-11, 1000}, {39, 500}, {39, 1000}, {1, 50}, {29, 500}},
                    {{12, 125}, {-79, 500}, {49, 500}, {-39, 1000}, {-11, 125}},
                    {{-17, 200}, {-169, 1000}, {153, 1000}, {-9, 50}, {21, 200}},
                    {{-17, 1000}, {27, 100}, {127, 1000}, {53, 500}, {107, 1000}},
                    {{39, 500}, {1, 250}, {-89, 500}, {67, 500}, {39, 1000}},
                    {{-199, 500}, {63, 1000}, {19, 1000}, {151, 1000}, {-61, 200}},
                    {{3, 40}, {9, 50}, {127, 1000}, {-79, 500}, {-141, 1000}},
                    {{-23, 250}, {-71, 1000}, {-8, 125}, {-97, 500}, {353, 1000}},
                    {{-12, 125}, {-101, 500}, {21, 1000}, {11, 100}, {-193, 1000}},
                    {{6, 125}, {-151, 1000}, {-37, 250}, {-1, 250}, {157, 1000}},
                    {{-7, 125}, {3, 25}, {3, 50}, {-29, 200}, {209, 1000}},
                    {{13, 40}, {-21, 200}, {-1, 1000}, {103, 1000}, {-7, 50}},
                    {{-211, 1000}, {-21, 100}, {73, 500}, {7, 500}, {-11, 125}},
                    {{-47, 500}, {9, 125}, {401, 1000}, {91, 1000}, {101, 1000}},
                    {{-3, 200}, {51, 500}, {-123, 500}, {3, 100}, {-24, 125}},
                    {{7, 50}, {-43, 1000}, {47, 200}, {167, 1000}, {63, 1000}},
                    {{-27, 250}, {-3, 40}, {281, 1000}, {61, 250}, {-37, 250}},
                    {{-31, 250}, {-1, 25}, {-59, 1000}, {31, 1000}, {129, 500}},
                    {{39, 250}, {31, 500}, {-181, 1000}, {39, 250}, {-143, 1000}},
                    {{0, 1}, {-9, 200}, {0, 1}, {53, 1000}, {43, 500}},
                    {{-41, 1000}, {1, 1000}, {-57, 1000}, {-163, 1000}, {-31, 200}},
                    {{43, 250}, {-1, 200}, {-13, 200}, {3, 500}, {-33, 500}},
                    {{-47, 1000}, {11, 200}, {-13, 200}, {11, 200}, {-293, 1000}},
                    {{-89, 1000}, {11, 125}, {1, 25}, {219, 1000}, {-37, 500}},
                    {{-1, 100}, {3, 100}, {-7, 100}, {163, 1000}, {27, 1000}}};

            path = generic_path<5>(reinterpret_cast<int32_t*>(vals), 50);
        }
    };

#include "accuracy_suite.ins"
#include "double_path_suite.ins"
#include "float_path_suite.ins"
#include "rational_path_suite.ins"
}

SUITE(brownian_path_10_2_10_tests)
{

    static const float expected_float_error = 2e-5f;
    static const double expected_double_error = 2e-13;

    //typedef GenericFixture<10, 2, 10> Fixture;

    struct Fixture : GenericFixture<10, 2, 10> {
        typedef GenericFixture<10, 2, 10> base;

        using base::rational_dense_framework;
        using base::rational_field;
        using base::rational_hybrid_framework;
        using base::rational_sparse_framework;

        using base::double_dense_framework;
        using base::double_field;
        using base::double_hybrid_framework;
        using base::double_sparse_framework;

        using base::float_dense_framework;
        using base::float_field;
        using base::float_hybrid_framework;
        using base::float_sparse_framework;

        Fixture()
        {
            int32_t vals[][10][2] = {
                    {{-269, 500}, {17, 1000}, {163, 200}, {53, 250}, {903, 1000}, {-63, 1000}, {561, 1000}, {-323, 1000}, {-117, 500}, {56, 125}},
                    {{-17, 200}, {-103, 1000}, {551, 1000}, {-11, 50}, {371, 1000}, {359, 1000}, {191, 500}, {-121, 500}, {61, 500}, {-63, 500}},
                    {{141, 250}, {-103, 1000}, {32, 125}, {521, 1000}, {37, 250}, {271, 1000}, {1, 100}, {49, 125}, {-187, 500}, {-111, 500}},
                    {{-451, 1000}, {-413, 1000}, {-453, 1000}, {-283, 500}, {51, 500}, {-27, 100}, {-451, 1000}, {131, 250}, {601, 1000}, {-213, 500}},
                    {{57, 1000}, {-141, 500}, {-137, 500}, {79, 250}, {171, 500}, {9, 500}, {209, 1000}, {107, 250}, {59, 250}, {13, 50}},
                    {{2, 5}, {-697, 1000}, {419, 1000}, {3, 5}, {-1, 40}, {21, 100}, {-9, 250}, {373, 1000}, {-23, 250}, {89, 250}},
                    {{51, 250}, {2, 25}, {1, 20}, {-9, 500}, {-97, 250}, {223, 1000}, {-269, 1000}, {39, 200}, {-31, 200}, {367, 1000}},
                    {{-517, 1000}, {-47, 1000}, {-69, 500}, {-37, 500}, {-289, 1000}, {-13, 200}, {-87, 1000}, {73, 200}, {63, 200}, {-21, 40}},
                    {{-303, 500}, {-201, 1000}, {179, 500}, {119, 250}, {-257, 1000}, {-123, 1000}, {-13, 1000}, {11, 250}, {-33, 500}, {-19, 250}},
                    {{-49, 1000}, {9, 100}, {27, 500}, {279, 1000}, {-57, 125}, {409, 1000}, {-13, 200}, {-29, 1000}, {41, 200}, {43, 100}}};

            path = generic_path<10>(reinterpret_cast<int32_t*>(vals), 10);
        }
    };

#include "accuracy_suite.ins"
#include "double_path_suite.ins"
#include "float_path_suite.ins"
#include "rational_path_suite.ins"
}
