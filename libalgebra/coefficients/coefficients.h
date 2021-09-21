//
// Created by sam on 25/06/2021.
//

#ifndef LIBALGEBRA_COEFFICIENTS_H
#define LIBALGEBRA_COEFFICIENTS_H

namespace alg {
namespace coefficients {

#define LIBALGEBRA_FIELD_GENERATE_BINARY(NAME, OP, RET_T, LHS_T, RHS_T)        \
  static LA_CONSTEXPR RET_T NAME(LHS_T const &lhs, RHS_T const &rhs) { return lhs OP rhs; } \
                                                                               \
  static LA_CONSTEXPR LHS_T &NAME##_inplace(LHS_T &lhs, RHS_T const &rhs) {                 \
    return (lhs OP## = rhs);                                                   \
  }

template <typename Scalar, typename Rational = Scalar> struct coefficient_field
{

    typedef Scalar SCA;
    typedef Rational RAT;
    typedef SCA S;
    typedef RAT Q;

    static const SCA one;
    static const SCA zero;
    static const SCA mone;

    static LA_CONSTEXPR SCA uminus(SCA arg) { return -arg; }

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

} // namespace coefficients
} // namespace alg

#endif // LIBALGEBRA_COEFFICIENTS_H
