//
// Created by sam on 21/01/2022.
//

#ifndef LIBALGEBRA_LIE_INNER_PRODUCT_H
#define LIBALGEBRA_LIE_INNER_PRODUCT_H

#include <libalgebra/coefficients/coefficients.h>
#include <libalgebra/coefficients/rational_coefficients.h>
#include <libalgebra/implementation_types.h>
#include <libalgebra/lie.h>
#include <libalgebra/lie_basis.h>
#include <libalgebra/operators/multi_linear_operators.h>
#include <libalgebra/tensor.h>
#include <libalgebra/utils.h>
#include <libalgebra/vectors/vectors.h>

#include <map>
#include <mutex>
#include <utility>

namespace alg {
namespace operators {

namespace dlt {

struct basic_tensor_inner_product {

    template<typename Tensor>
    typename Tensor::SCA operator()(const Tensor& t1, const Tensor& t2) const
    {
        typename Tensor::SCA result;
        for (const auto kvp : t1) {
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
    using arg_type = alg::lie<Coeffs, Width, Depth, ArgVecType, ArgArgs...>;

    using cache_coeff_t = coefficients::rational_field;
    using cache_vec_type = alg::lie<cache_coeff_t, Width, Depth, vectors::sparse_vector>;
    using tensor_type = alg::free_tensor<cache_coeff_t, Width, Depth, vectors::sparse_vector>;

    static std::map<key_pair, cache_vec_type> cache;

    typename cache_coeff_t::S key_level_impl(const key_type& key1, const key_type& key2) const
    {
        maps<cache_coeff_t, Width, Depth, cache_vec_type, tensor_type> maps;

        const auto t1 = maps.l2t(cache_vec_type(key1, cache_coeff_t::one));
        const auto t2 = maps.l2t(cache_vec_type(key2, cache_coeff_t::one));

        return TensorInnerProduct(t1, t2);
    }

    const typename cache_coeff_t::S& key_level(const key_type& key1, const key_type& key2) const
    {
        basis_type& basis = cache_vec_type::basis;

        key_pair p(key1, key2);

        const auto cache_it = cache.find(p);
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
        for (const auto& a1_it : arg1) {
            for (const auto& a2_it : arg2) {
                result += scal(key_level(a1_it.key(), a2_it.key())) * a1_it.value() * a2_it.value();
            }
        }
        return result;
    }
};

}// namespace dlt

}// namespace operators
}// namespace alg

#endif//LIBALGEBRA_LIE_INNER_PRODUCT_H
