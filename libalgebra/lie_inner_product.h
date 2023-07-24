//
// Created by sam on 21/01/2022.
//

#ifndef LIBALGEBRA_LIE_INNER_PRODUCT_H
#define LIBALGEBRA_LIE_INNER_PRODUCT_H

#include "implementation_types.h"

#include <map>
#include <mutex>
#include <utility>

#include "vectors.h"
#include "rational_coefficients.h"
#include "lie.h"
#include "tensor.h"
#include "utils.h"


namespace alg {
namespace operators {

namespace dtl {

struct basic_tensor_inner_product {

    template<typename Tensor>
    typename Tensor::SCA operator()(const Tensor& t1, const Tensor& t2) const
    {
        typename Tensor::SCA result;
        for (auto kvp : t1) {
            result += kvp.value() * t2[kvp.key()];
        }
        return result;
    }
};

template<
        DEG Width,
        DEG Depth,
        typename TensorInnerProduct = basic_tensor_inner_product>
class lie_inner_product_impl
{
    using basis_type = lie_basis<Width, Depth>;
    using key_type = typename basis_type::KEY;
    using key_pair = std::pair<key_type, key_type>;

    template<typename Coeffs, template<typename, typename, typename...> class ArgVecType, typename... ArgArgs>
    using arg_type = lie<Coeffs, Width, Depth, ArgVecType, ArgArgs...>;

    using cache_coeff_t = coefficients::rational_field;
    using cache_vec_type = lie<cache_coeff_t, Width, Depth, vectors::sparse_vector>;
    using tensor_type = free_tensor<cache_coeff_t, Width, Depth, vectors::sparse_vector>;

    TensorInnerProduct tensor_ip;

    static std::map<key_pair, typename cache_coeff_t::S> cache;

    typename cache_coeff_t::S key_level_impl(const key_type& key1, const key_type& key2) const
    {
        maps<cache_coeff_t, Width, Depth, tensor_type, cache_vec_type> maps;

        const auto t1 = maps.l2t(cache_vec_type(key1, cache_coeff_t::one));
        const auto t2 = maps.l2t(cache_vec_type(key2, cache_coeff_t::one));

        return tensor_ip(t1, t2);
    }

    const typename cache_coeff_t::S& key_level(const key_type& key1, const key_type& key2) const
    {
//        basis_type& basis = cache_vec_type::basis;

        key_pair p(key1, key2);

        auto cache_it = cache.find(p);
        if (cache_it != cache.end()) {
            return cache_it->second;
        }

        return cache[p] = key_level_impl(key1, key2);
    }

public:
    template<
            typename Coeffs,
            template<typename, typename, typename...> class Arg1VecType,
            typename... Arg1Args,
            template<typename, typename, typename...> class Arg2VecType,
            typename... Arg2Args>
    typename Coeffs::S operator()(
            const arg_type<Coeffs, Arg1VecType, Arg1Args...>& arg1,
            const arg_type<Coeffs, Arg2VecType, Arg2Args...>& arg2) const
    {
        using scal = typename Coeffs::S;
        scal result;
        for (auto& a1_it : arg1) {
            for (auto& a2_it : arg2) {
                result += scal(key_level(a1_it.key(), a2_it.key())) * a1_it.value() * a2_it.value();
            }
        }
        return result;
    }
};

template<
        DEG Width,
        DEG Depth,
        typename TensorInnerProduct>
std::map<std::pair<typename lie_basis<Width, Depth>::KEY, typename lie_basis<Width, Depth>::KEY>,
         typename coefficients::rational_field::S>
        lie_inner_product_impl<Width, Depth, TensorInnerProduct>::cache;

}// namespace dtl

template<
        typename Coeffs,
        DEG Width,
        DEG Depth,
        template<typename, typename, typename...> class ArgVecType,
        typename... ArgArgs>
using lie_inner_product = multi_linear_operator<
        dtl::lie_inner_product_impl<Width, Depth>,
        lie<Coeffs, Width, Depth, ArgVecType, ArgArgs...>,
        lie<Coeffs, Width, Depth, ArgVecType, ArgArgs...>>;

}// namespace operators
}// namespace alg

#endif//LIBALGEBRA_LIE_INNER_PRODUCT_H
