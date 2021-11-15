//
// Created by sam on 15/10/2021.
//

#ifndef LIBALGEBRA_RANDOM_VECTOR_GENERATOR_H
#define LIBALGEBRA_RANDOM_VECTOR_GENERATOR_H

#include <boost/random.hpp>
#include <utility>
#include <vector>

namespace la_testing {

template<typename Vector, typename CoeffDist = boost::random::uniform_real_distribution<typename Vector::SCALAR>, typename SkipDist = void>
struct random_vector_generator {
    using basis_type = typename Vector::BASIS;
    using key_type = typename basis_type::KEY;
    using scalar_type = typename Vector::SCALAR;
    using kv_pair = std::pair<key_type, scalar_type>;
    using kv_vec = std::vector<kv_pair>;
    using size_type = typename SkipDist::result_type;

    template<typename... Args>
    explicit random_vector_generator(const SkipDist& skip, Args&&... args)
        : m_coeff_dist(std::forward<Args>(args)...), m_skip_dist(skip)
    {}

    template<typename Rng>
    Vector operator()(Rng& rng)
    {
        const basis_type& basis = Vector::basis;
        Vector tmp;

        size_type skip = m_skip_dist(rng);
        for (key_type k = basis.begin(); k != basis.end(); k = basis.nextkey(k)) {
            if (skip == 0) {
                tmp.add_scal_prod(k, m_coeff_dist(rng));
                skip = m_skip_dist(rng);
            }
            else {
                --skip;
            }
        }

        return tmp;
    }

private:
    CoeffDist m_coeff_dist;
    SkipDist m_skip_dist;
};

/**
 * @brief Generate random vectors in a dense way
 * @tparam Vector Vector type
 * @tparam CoeffDist Distribution function of coefficients
 */
template<typename Vector, typename CoeffDist>
struct random_vector_generator<Vector, CoeffDist, void> {
    using basis_type = typename Vector::BASIS;
    using key_type = typename basis_type::KEY;
    using scalar_type = typename Vector::SCALAR;
    using kv_pair = std::pair<key_type, scalar_type>;
    using kv_vec = std::vector<kv_pair>;

    template<typename... Args>
    explicit random_vector_generator(Args&&... args) : m_coeff_dist(std::forward<Args>(args)...)
    {}

    template<typename Rng>
    Vector operator()(Rng& rng)
    {
        const basis_type& basis = Vector::basis;
        Vector tmp;

        for (key_type k = basis.begin(); k != basis.end(); k = basis.nextkey(k)) {
            tmp.add_scal_prod(k, m_coeff_dist(rng));
        }

        return tmp;
    }

private:
    CoeffDist m_coeff_dist;
};

}// namespace la_testing

#endif//LIBALGEBRA_RANDOM_VECTOR_GENERATOR_H

