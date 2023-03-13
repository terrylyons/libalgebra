//
// Created by sam on 12/02/2021.
//

#include <libalgebra/algebra.h>
#include "simple_basis.h"

template <typename Basis>
class pointwise_multiplier : alg::multiplier_base<pointwise_multiplier<Basis>, Basis>
{
    using base = alg::multiplier_base<pointwise_multiplier, Basis>;
    friend base;

public:

    using basis_type = Basis;
    using key_type = typename Basis::KEY;
    using typename base::pair_type;
    using result_type = typename base::inner_result_type;
    using typename base::argument_type;


    result_type operator()(argument_type lhs, argument_type rhs) const
    {
        if (lhs != rhs) {
            return {};
        }
        return {{lhs, 1}};
    }

};

template <typename Basis>
using pointwise_multiplication = alg::base_multiplication<pointwise_multiplier<Basis>>;
