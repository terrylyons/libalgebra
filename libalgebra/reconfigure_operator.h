//
// Created by user on 16/03/23.
//

#ifndef LIBALGEBRA_LIBALGEBRA_RECONFIGURE_OPERATOR_H_
#define LIBALGEBRA_LIBALGEBRA_RECONFIGURE_OPERATOR_H_

#include "algebra.h"
#include "implementation_types.h"
#include "lie.h"
#include "operators.h"
#include "tensor.h"

#include <type_traits>
#include <utility>

namespace alg {
namespace operators {

template<DEG NewDepth>
struct change_depth {
    static constexpr DEG new_depth = NewDepth;
};

template<typename NewCoeffs>
struct change_coeffs {
    using new_coeffs = NewCoeffs;
};

template<template<typename, typename, typename...> class NewVType, typename... NewArgs>
struct change_storage {
    template<typename B, typename C, typename... Args>
    using new_storage = NewVType<B, C, Args...>;
};

namespace dtl {


template<typename R, typename Base>
class recalibrate;

template<DEG NewDepth, typename Base>
class recalibrate<change_depth<NewDepth>, Base> : Base
{
    using old_type = typename Base::type;

    template<typename C, DEG W, DEG OD, template<typename, typename, typename...> class V, template<DEG, DEG> class M, typename... Args>
    static free_tensor<C, W, NewDepth, V, M, Args...> f(free_tensor<C, W, OD, V, M, Args...>& a);

    template<typename C, DEG W, DEG OD, typename... Args>
    static shuffle_tensor<C, W, NewDepth, Args...> f(shuffle_tensor<C, W, OD, Args...>& a);

    template<typename C, DEG W, DEG OD, template<typename, typename, typename...> class V, typename... Args>
    static lie<C, W, NewDepth, V, Args...> f(lie<C, W, OD, V, Args...>& a);

public:
    using type = decltype(f(std::declval<old_type&>()));
};

template<typename NewCoeffs, typename Base>
class recalibrate<change_coeffs<NewCoeffs>, Base> : Base
{
    using old_type = typename Base::type;

    template<typename C, DEG W, DEG D, template<typename, typename, typename...> class V, template<DEG, DEG> class M, typename... Args>
    static free_tensor<NewCoeffs, W, D, V, M, Args...> f(free_tensor<C, W, D, V, M, Args...>& a);

    template<typename C, DEG W, DEG D, typename... Args>
    static shuffle_tensor<NewCoeffs, W, D, Args...> f(shuffle_tensor<C, W, D, Args...>& a);

    template<typename C, DEG W, DEG D, template<typename, typename, typename...> class V, typename... Args>
    static lie<NewCoeffs, W, D, V, Args...> f(lie<C, W, D, V, Args...>& a);


public:
    using type = decltype(f(std::declval<old_type&>()));
};

template<template <typename, typename, typename...> class NewVType, typename Base, typename... NewVArgs>
class recalibrate<change_storage<NewVType, NewVArgs...>, Base> : Base
{
    using old_type = typename Base::type;

    template<typename C, DEG W, DEG D, template<typename, typename, typename...> class V, template<DEG, DEG> class M, typename... Args>
    static free_tensor<C, W, D, NewVType, M, NewVArgs...> f(free_tensor<C, W, D, V, M, Args...>& a);

    template<typename C, DEG W, DEG D, typename... Args>
    static shuffle_tensor<C, W, D, Args...> f(shuffle_tensor<C, W, D, Args...>& a);

    template<typename C, DEG W, DEG D, template<typename, typename, typename...> class V, typename... Args>
    static lie<C, W, D, NewVType, NewVArgs...> f(lie<C, W, D, V, Args...>& a);


public:
    using type = decltype(f(std::declval<old_type&>()));
};


template <typename Algebra, typename... Changes>
struct recalibrate_chain;

template <typename Algebra, typename FirstChange, typename... Changes>
struct recalibrate_chain<Algebra, FirstChange, Changes...>
    : recalibrate<FirstChange, recalibrate_chain<Algebra, Changes...>>
{};

template <typename Algebra>
struct recalibrate_chain<Algebra> {
    using type = Algebra;
};

template <typename InputType, typename... Changes>
using recalibrated_t = typename recalibrate_chain<InputType, Changes...>::type;


template <typename Algebra, typename... Changes>
class recalibrate_operator_impl {

    template <template <DEG, DEG> class B, DEG W, DEG D1,
             DEG D2,
             typename C1,
             typename C2>
    void construct_impl(vectors::dense_vector<B<W, D1>, C1>& dst,
                        const vectors::dense_vector<B<W, D2>, C2>& src) {
        auto max_deg = std::min(D1, src.degree());
        dst.resize_to_degree(max_deg);

        auto* dst_ptr = dst.as_mut_ptr();
        const auto* src_ptr = src.as_ptr();

        auto dim = dst.dimension();

        for (DIMN i=0; i<dim; ++i) {
            dst_ptr[i] = static_cast<typename C2::S>(src_ptr[i]);
        }

    }

    template <template <DEG, DEG> class B, DEG W, DEG D1,
             DEG D2,
             typename C1,
             typename C2,
             template <typename, typename, typename...> class V1,
             template <typename, typename, typename...> class V2>
    void construct_impl(V1<B<W, D1>, C1>& dst,
                        const V2<B<W, D2>, C2>& src) {
        for (auto&& it : src) {
            auto key = it.key();
            if (src.degree(key) <= D1) {
                dst[it.key()] = static_cast<typename C1::S>(src.value());
            }
        }
    }




public:

    using argument_type = Algebra;
    using result_type = recalibrated_t<Algebra, Changes...>;



    result_type operator()(const argument_type& arg) const {
        result_type result;
        construct_impl(result.base_vector(), arg.base_vector());
        return result;
    }


};




}// namespace dtl

template<typename Algebra, typename... Options>
using recalibrate_operator = linear_operator<
        dtl::recalibrate_operator_impl<Algebra, Options...>,
        Algebra, dtl::recalibrated_t<Algebra, Options...>>;

}// namespace operators
}// namespace alg

#endif//LIBALGEBRA_LIBALGEBRA_RECONFIGURE_OPERATOR_H_
