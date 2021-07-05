//
// Created by sam on 25/06/2021.
//

#ifndef LIBALGEBRA_COEFFICIENTS_H
#define LIBALGEBRA_COEFFICIENTS_H

namespace alg {
namespace coefficients {

#define LIBALGEBRA_FIELD_GENERATE_BINARY(NAME, OP, RET_T, LHS_T, RHS_T) \
    static RET_T NAME(LHS_T lhs, RHS_T rhs) {                           \
        return lhs OP rhs;                                              \
    }

template <typename Scalar, typename Rational=Scalar>
struct coefficient_field {

    typedef Scalar SCA;
    typedef Rational RAT;

    static const SCA one;
    static const SCA zero;
    static const SCA mone;

    static SCA uminus(SCA arg) {
        return -arg;
    }

    LIBALGEBRA_FIELD_GENERATE_BINARY(add, +, SCA, SCA, SCA);

    LIBALGEBRA_FIELD_GENERATE_BINARY(sub, -, SCA, SCA, SCA);

    LIBALGEBRA_FIELD_GENERATE_BINARY(mul, *, SCA, SCA, SCA);

    LIBALGEBRA_FIELD_GENERATE_BINARY(div, /, SCA, SCA, RAT);


};

template <typename Scalar, typename Rational> const Scalar coefficient_field<Scalar, Rational>::one(1);

template <typename Scalar, typename Rational> const Scalar coefficient_field<Scalar, Rational>::zero(0);

template <typename Scalar, typename Rational> const Scalar coefficient_field<Scalar, Rational>::mone(-1);


#undef LIBALGEBRA_FIELD_GENERATE_BINARY


typedef coefficient_field<double> double_field;
typedef coefficient_field<float> float_field;

} // namespace coeffcients
} // namespace alg


#endif //LIBALGEBRA_COEFFICIENTS_H
