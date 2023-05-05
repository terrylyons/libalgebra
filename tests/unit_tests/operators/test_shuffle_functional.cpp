//
// Created by sam on 01/11/2021.
//

#include <UnitTest++.h>

#include <libalgebra/alg_types.h>
#include <libalgebra/functionals.h>
#include <libalgebra/operators.h>

#include "../../common/random_vector_generator.h"
#include "../../common/rng.h"
#include <memory>

using namespace alg;

template<typename ShuffleTensor>
class shuffle_tensor_implementation
{
    using basis_t = typename ShuffleTensor::BASIS;
    using ft_basis_t = free_tensor_basis<basis_t::s_no_letters, basis_t::s_max_degree>;
    using coeff_t = typename ShuffleTensor::coefficient_field;
    using scalar_t = typename coeff_t::S;

public:
    explicit shuffle_tensor_implementation(ShuffleTensor&& t) : m_data(std::move(t))
    {}

    template<typename... Args>
    explicit shuffle_tensor_implementation(Args&&... args) : m_data(std::forward<Args>(args)...)
    {}

    template<typename ArgumentType>
    scalar_t operator()(const ArgumentType& arg) const
    {
        scalar_t result(coeff_t::zero);
        for (auto& it : m_data) {
            coeff_t::add_inplace(result, coeff_t::mul(it.value(), arg[it.key()]));
        }
        return result;
    }

private:
    ShuffleTensor m_data;
};

SUITE(shuffle_tensor_functional)
{

    struct fixture : alg_types<5, 5, DPReal> {
        using functional_t = operators::linear_functional<shuffle_tensor_implementation<SHUFFLE_TENSOR>>;
    };

    TEST_FIXTURE(fixture, test_bracket_zeros)
    {
        TENSOR ftensor;
        functional_t func;

        CHECK_EQUAL(0.0, func(ftensor));
    }

    TEST_FIXTURE(fixture, test_bracket_identity)
    {
        TENSOR ftensor(S(1));
        functional_t func(S(1));

        CHECK_EQUAL(1.0, func(ftensor));
    }

    TEST_FIXTURE(fixture, test_bracket_deg1)
    {
        std::vector<double> fdata{0.0, 1.0, 2.0, 3.0};
        std::vector<double> sdata{0.0, 1.0, 1.0, 1.0};

        TENSOR ftensor(fdata.data(), fdata.data() + fdata.size());
        functional_t func(sdata.data(), sdata.data() + sdata.size());

        CHECK_EQUAL(6.0, func(ftensor));
    }

    TEST_FIXTURE(fixture, test_bracket_deg1_mismatch)
    {
        std::vector<double> fdata{0.0, 1.0, 2.0, 3.0};
        std::vector<double> sdata{0.0, 0.0, 1.0, 0.0};

        TENSOR ftensor(fdata.data(), fdata.data() + fdata.size());
        functional_t func(sdata.data(), sdata.data() + sdata.size());

        CHECK_EQUAL(2.0, func(ftensor));
    }

    TEST_FIXTURE(fixture, test_bracket_deg1_complete_mismatch)
    {
        std::vector<double> fdata{0.0, 1.0, 2.0, 3.0};
        std::vector<double> sdata{1.0, 0.0, 0.0, 0.0};

        TENSOR ftensor(fdata.data(), fdata.data() + fdata.size());
        functional_t func(sdata.data(), sdata.data() + sdata.size());

        CHECK_EQUAL(0.0, func(ftensor));
    }
}

SUITE(shuffle_tensor_library_implementation)
{

    struct fixture : alg_types<5, 5, DPReal> {
        using functional_t = operators::shuffle_tensor_functional<SHUFFLE_TENSOR, TENSOR>;
    };

    TEST_FIXTURE(fixture, test_bracket_zeros)
    {
        TENSOR ftensor;
        functional_t func;

        CHECK_EQUAL(0.0, func(ftensor));
    }

    TEST_FIXTURE(fixture, test_bracket_identity)
    {
        TENSOR ftensor(S(1));
        functional_t func(S(1));

        CHECK_EQUAL(1.0, func(ftensor));
    }

    TEST_FIXTURE(fixture, test_bracket_deg1)
    {
        std::vector<double> fdata{0.0, 1.0, 2.0, 3.0};
        std::vector<double> sdata{0.0, 1.0, 1.0, 1.0};

        TENSOR ftensor(fdata.data(), fdata.data() + fdata.size());
        functional_t func(sdata.data(), sdata.data() + sdata.size());

        CHECK_EQUAL(6.0, func(ftensor));
    }

    TEST_FIXTURE(fixture, test_bracket_deg1_mismatch)
    {
        std::vector<double> fdata{0.0, 1.0, 2.0, 3.0};
        std::vector<double> sdata{0.0, 0.0, 1.0, 0.0};

        TENSOR ftensor(fdata.data(), fdata.data() + fdata.size());
        functional_t func(sdata.data(), sdata.data() + sdata.size());

        CHECK_EQUAL(2.0, func(ftensor));
    }

    TEST_FIXTURE(fixture, test_bracket_deg1_complete_mismatch)
    {
        std::vector<double> fdata{0.0, 1.0, 2.0, 3.0};
        std::vector<double> sdata{1.0, 0.0, 0.0, 0.0};

        TENSOR ftensor(fdata.data(), fdata.data() + fdata.size());
        functional_t func(sdata.data(), sdata.data() + sdata.size());

        CHECK_EQUAL(0.0, func(ftensor));
    }


}
