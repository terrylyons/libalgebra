//
// Created by sam on 24/01/2022.
//

#ifndef LIBALGEBRA_RANDOM_COEFFS_H
#define LIBALGEBRA_RANDOM_COEFFS_H

#include <libalgebra/libalgebra.h>
#include <libalgebra/coefficients.h>
#include <libalgebra/detail/meta.h>

#include <boost/random.hpp>
#include <random>

namespace la_testing {

template<typename Coeffs, template<typename> class Dist = boost::random::uniform_real_distribution>
struct random_coeff_generator {
    using base_coeff = typename alg::utils::scalar_base<typename Coeffs::SCA>::type;
    using distribution = Dist<base_coeff>;

    distribution base_dist;

    template<typename... Args>
    explicit random_coeff_generator(Args&&... args) : base_dist(std::forward<Args>(args)...)
    {}

    template<typename U, typename RNG>
    typename std::enable_if<std::is_same<base_coeff, U>::value, base_coeff>::type
    generate(RNG& rng)
    {
        return base_dist(rng);
    }

    template<typename U, typename RNG>
    typename std::enable_if<alg::utils::is_vector_type<U>::value, U>::type
    generate(RNG& rng)
    {
        using vector_type = U;
        auto& basis = vector_type::basis;

        vector_type result;
        for (auto k = basis.begin(); k != basis.end(); k = basis.nextkey(k)) {
            result.add_scal_prod(k, base_dist(rng));
        }

        return result;
    }

    template<typename RNG>
    auto operator()(RNG& rng) -> decltype(generate<typename Coeffs::SCA>(rng))
    {
        return generate<typename Coeffs::SCA>(rng);
    }
};

}// namespace la_testing

#endif//LIBALGEBRA_RANDOM_COEFFS_H
