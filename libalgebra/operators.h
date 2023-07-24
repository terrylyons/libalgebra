//
// Created by sam on 29/10/2021.
//

#ifndef LIBALGEBRA_OPERATORS_H
#define LIBALGEBRA_OPERATORS_H

#include <type_traits>
#include <utility>

#include "composition_operator.h"
#include "scalar_multiply_operator.h"
#include "sum_operator.h"
#include "vector_bundle.h"



namespace alg {
namespace operators {

namespace dtl {

template <typename Impl, typename ArgBase, typename ArgFibre>
struct result_is_not_scalar_bundle {
    using base_result_t = decltype(std::declval<Impl>()(std::declval<const ArgBase&>()));
    using fibre_result_t = decltype(std::declval<Impl>()(std::declval<const ArgFibre&>()));


    static constexpr bool value = !std::is_same<base_result_t, typename ArgBase::SCALAR>::value
            && !std::is_same<fibre_result_t, typename ArgFibre::SCALAR>::value;

};

}

template<typename Impl>
class linear_operator : protected Impl
{

public:

protected:
    using implementation_type = Impl;

public:
    using implementation_type::implementation_type;

    using implementation_type::operator();

//    template <typename ArgBase, typename ArgFibre,
//             typename=std::enable_if_t<dtl::result_is_not_scalar_bundle<Impl, ArgBase, ArgFibre>::value>>
    template <typename ArgBase, typename ArgFibre>
    auto operator()(const vector_bundle<ArgBase, ArgFibre>& arg)
            -> vector_bundle<decltype(Impl::operator()(arg.base())),
                    decltype(Impl::operator()(arg.fibre()))>
    {
        return {Impl::operator()(arg.base()), Impl::operator()(arg.fibre())};
    }

protected:
    explicit linear_operator(Impl&& impl) : Impl(std::forward<Impl>(impl))
    {}

public:
    template<typename OtherImpl>
    linear_operator<sum_operator<Impl, OtherImpl>>
    operator+(const linear_operator<OtherImpl>& other) const
    {
        return {*this, other};
    }

    template<typename OtherImpl, typename NewArgumentType>
    linear_operator<composition_operator<Impl, OtherImpl>> friend
    operator*(const linear_operator& outer,
              const linear_operator<OtherImpl>& inner)
    {
        using composition_type = composition_operator<Impl, OtherImpl>;
        using return_type = linear_operator<composition_type>;

        return return_type(composition_type(static_cast<const OtherImpl&>(inner)), static_cast<const Impl&>(outer));
    }
};

namespace dtl {

template<typename Algebra>
class left_multiplication_operator_impl
{
protected:
    using shuffle_algebra_t = Algebra;

public:
    template<typename... Args>
    explicit left_multiplication_operator_impl(Args&&... args)
        : m_lhs(std::forward<Args>(args)...)
    {}

    explicit left_multiplication_operator_impl(shuffle_algebra_t&& alg) : m_lhs(alg)
    {}

    shuffle_algebra_t operator()(const shuffle_algebra_t& arg) const
    {
        return m_lhs * arg;
    }

private:
    shuffle_algebra_t m_lhs;
};

template<typename Algebra>
class right_multiplication_operator_impl
{
protected:
    using shuffle_algebra_t = Algebra;

public:
    template<typename... Args>
    explicit right_multiplication_operator_impl(Args&&... args)
        : m_rhs(std::forward<Args>(args)...)
    {}

    explicit right_multiplication_operator_impl(shuffle_algebra_t&& alg) : m_rhs(alg)
    {}

    shuffle_algebra_t operator()(const shuffle_algebra_t& arg) const
    {
        return arg * m_rhs;
    }

private:
    shuffle_algebra_t m_rhs;
};

}//namespace dtl

template<typename Algebra>
using left_multiplication_operator = linear_operator<
        dtl::left_multiplication_operator_impl<Algebra>>;

template<typename Algebra>
using right_multiplication_operator = linear_operator<
        dtl::right_multiplication_operator_impl<Algebra>>;

namespace dtl {

template<typename FreeTensor>
class adjoint_of_left_multiplication_operator_impl
{
protected:
    using tensor_algebra_t =  FreeTensor;

public:
    template<typename... Args>
    explicit adjoint_of_left_multiplication_operator_impl(Args&&... args)
        : m_lhs(std::forward<Args>(args)...)
    {}

    explicit adjoint_of_left_multiplication_operator_impl(tensor_algebra_t&& alg) : m_lhs(alg)
    {}

    template <typename ShuffleTensor>
    std::enable_if_t<ShuffleTensor::BASIS::s_no_letters == FreeTensor::BASIS::s_no_letters, ShuffleTensor>
    operator()(const ShuffleTensor& arg) const
    {
        ShuffleTensor result;
        eval(result.base_vector(), arg.base_vector(), m_lhs.base_vector());
        return result;
    }

private:
    tensor_algebra_t m_lhs;

    template<typename Arg>
    static Arg shift_down(const Arg& sh, typename Arg::KEY word)
    {
        Arg result(sh);
        Arg working;
        while (word.size()) {
            auto letter = word.lparent();
            word = word.rparent();
            for (auto pr : result) {
                if (pr.key().size() > 0 && pr.key().lparent() == letter)
                    working[pr.key().rparent()] = result[pr.key()];
            }
            result.swap(working);
            working.clear();
        }
        return result;
    }

    template<typename Arg, typename FTensor>
    static void eval(Arg& result, const Arg& arg, const FTensor& param)
    {
        using s_t = typename Arg::SCALAR;
        for (auto&& pr : param) {
            const auto& val = pr.value();
            Arg::template apply_inplace_flat_binary_op(result, shift_down(arg, pr.key()),
                                                       [&val](const s_t& l, const s_t& r) { return l + val * r; });
        }
    }

    template<typename S, typename T>
    static void eval_single(S* LA_RESTRICT out,
                            const S* LA_RESTRICT in,
                            const T* LA_RESTRICT param,
                            const DIMN* powers,
                            const DIMN* sizes,
                            DEG param_deg,
                            DEG arg_deg)
    {
        if (param_deg < arg_deg) {
            auto* out_ptr = out;
            for (DEG d = param_deg; d <= arg_deg; ++d) {
                auto result_deg = d - param_deg;

                for (DIMN param_idx = 0; param_idx < powers[param_deg]; ++param_idx) {
                    const auto* in_ptr = in + param_idx * powers[result_deg];
                    const auto& param_val = param[param_idx];

                    for (DIMN i = 0; i < powers[result_deg]; ++i) {
                        out_ptr[i] += param_val * in_ptr[i];
                    }
                }

                in += powers[d];
                out_ptr += powers[result_deg];
            }
        }
        else {
            auto& unit = out[0];

            for (DIMN i=0; i<powers[param_deg]; ++i) {
                unit += param[i]*in[i];
            }

        }
    }

    template <typename ShB, typename C, typename FTensor>
    static void eval(vectors::dense_vector<ShB, C>& result, const vectors::dense_vector<ShB, C>& arg, const FTensor& param) {
        using tsi = typename ShB::SIZE_INFO;

        const auto& basis = result.basis;

        auto arg_deg = arg.degree();
        auto param_deg = param.degree();
        result.resize_to_degree(arg_deg);

        for (auto&& item : param) {
            const auto prefix_deg = item.first.size();
            auto idx = basis.key_to_index(item.first);
            idx -= basis.start_of_degree(prefix_deg);

            for (DEG deg = prefix_deg; deg <= arg_deg; ++deg) {
                const auto suffix_deg = deg - prefix_deg;
                auto* optr = arg.as_mut_ptr() + basis.start_of_degree(suffix_deg);
                const auto* iptr = arg.as_ptr() + idx * tsi::powers[suffix_deg];

                for (DIMN i=0; i<tsi::powers[suffix_deg]; ++i) {
                    optr[i] += item.second * iptr[i];
                }
            }

        }


    }

    template<typename ShB, typename FB, typename C>
    static void eval(vectors::dense_vector<ShB, C>& result, const vectors::dense_vector<ShB, C>& arg, const vectors::dense_vector<FB, C>& param)
    {
        using tsi = typename ShB::SIZE_INFO;

        auto arg_deg = arg.degree();
        auto param_deg = param.degree();
        DEG lower_deg = std::min(arg_deg, param_deg);
        result.resize_to_degree(arg_deg);

        auto* out_ptr = result.as_mut_ptr();
        const auto* arg_ptr = arg.as_ptr();
        const auto* param_ptr = param.as_ptr();

        const auto& param_unit = param_ptr[0];
        if (param_unit != typename C::S(0)) {
            for (DIMN i = 0; i < arg.dimension(); ++i) {
                out_ptr[i] = param_unit * arg_ptr[i];
            }
        }

        for (DEG prefix_deg = 1; prefix_deg <= lower_deg; ++prefix_deg) {
            arg_ptr += tsi::powers[prefix_deg - 1];
            param_ptr += tsi::powers[prefix_deg - 1];

            eval_single(out_ptr, arg_ptr, param_ptr, tsi::powers.data(), tsi::degree_sizes.data(), prefix_deg, arg_deg);
        }
    }
};

}// namespace dtl

template<typename TensorAlgebra>
using adjoint_of_left_multiplication_operator = linear_operator<
        dtl::adjoint_of_left_multiplication_operator_impl<TensorAlgebra>>;

}// namespace operators
}// namespace alg

#include "free_extension.h"
#include "tensor_operator.h"

#endif//LIBALGEBRA_OPERATORS_H
