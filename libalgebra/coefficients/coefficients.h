//
// Created by sam on 25/06/2021.
//

#ifndef LIBALGEBRA_COEFFICIENTS_H
#define LIBALGEBRA_COEFFICIENTS_H

#include <libalgebra/implementation_types.h>
#include <libalgebra/utils/meta.h>

#include <type_traits>

namespace alg {
namespace coefficients {

#define LIBALGEBRA_FIELD_GENERATE_BINARY(NAME, OP, RET_T, LHS_T, RHS_T)                                     \
    static constexpr RET_T NAME(LHS_T const& lhs, RHS_T const& rhs) noexcept(noexcept(lhs OP rhs))          \
    {                                                                                                       \
        return lhs OP rhs;                                                                                  \
    }                                                                                                       \
                                                                                                            \
    static constexpr LHS_T& NAME##_inplace(LHS_T& lhs, RHS_T const& rhs) noexcept(noexcept(lhs OP## = rhs)) \
    {                                                                                                       \
        return (lhs OP## = rhs);                                                                            \
    }

template<typename Scalar, typename Rational>
struct coefficient_ring {

    static constexpr bool is_unital = std::is_constructible<Scalar, typename alg::utils::scalar_base<Scalar>::type>::value;

    typedef Scalar SCA;
    typedef Rational RAT;
    typedef SCA S;
    typedef RAT Q;

    static const SCA zero;

    static const typename std::enable_if<is_unital, SCA>::type one;
    static const typename std::enable_if<is_unital, SCA>::type mone;

    static constexpr Scalar from(DEG degree) noexcept(noexcept(Scalar(degree)))
    {
        return static_cast<Scalar>(degree);
    }

    template<typename T>
    static constexpr Scalar from(T arg) noexcept(noexcept(Scalar(arg)))
    {
        return static_cast<Scalar>(arg);
    }

    static constexpr SCA uminus(SCA arg) noexcept(noexcept(-arg))
    {
        return -arg;
    }

    LIBALGEBRA_FIELD_GENERATE_BINARY(add, +, SCA, SCA, SCA);

    LIBALGEBRA_FIELD_GENERATE_BINARY(sub, -, SCA, SCA, SCA);

    LIBALGEBRA_FIELD_GENERATE_BINARY(mul, *, SCA, SCA, SCA);

    LIBALGEBRA_FIELD_GENERATE_BINARY(div, /, SCA, SCA, RAT);
};

/*
 * The zero element should always exist and should be the result of default initialisation.
 */
template<typename Scalar, typename Rational>
const Scalar coefficient_ring<Scalar, Rational>::zero;

/*
 * One and minus one only exist when the ring is unital.
 */
template<typename Scalar, typename Rational>
const typename std::enable_if<
        coefficient_ring<Scalar, Rational>::is_unital,
        Scalar>::type coefficient_ring<Scalar, Rational>::one(1);

template<typename Scalar, typename Rational>
const typename std::enable_if<
        coefficient_ring<Scalar, Rational>::is_unital,
        Scalar>::type coefficient_ring<Scalar, Rational>::mone(-1);

/**
 * @brief A field (or more generally, ring) that can be used as coefficients in a vector
 *
 * A coefficient field type encapsulates the properties that are required to build
 * a valid vector space (or module if the coefficients do not actually form a field).
 *
 * The most important items are, of course, the declaration of the scalar type S
 * and rational type Q. These should be the same if the coefficients do form a field,
 * but can be different in general. The class also hold static members that store
 * the basic elements: one, minus one, and zero. Currently these are also stored by
 * the vector classes, but in the future the plan is to remove those static members
 * and use these values instead.
 *
 * Finally, the class contains some wrapper functions around the arithmetic operations
 * for the field. At the moment, these simply invoke the underlying operations but in
 * the future this mechanism can be changed to allow for different operations beyond
 * the built in arithmetic for type like floats and doubles.
 *
 * @tparam Scalar Type for scalar (coefficient) values
 */
template<typename Scalar>
struct coefficient_field : public coefficient_ring<Scalar, Scalar> {

    using ring_type = coefficient_ring<Scalar, Scalar>;

    static_assert(ring_type::is_unital, "a field must be unital");

    using typename ring_type::Q;
    using typename ring_type::RAT;
    using typename ring_type::S;
    using typename ring_type::SCA;


};

#undef LIBALGEBRA_FIELD_GENERATE_BINARY

typedef coefficient_field<double> double_field;
typedef coefficient_field<float> float_field;

}// namespace coefficients
}// namespace alg

#endif// LIBALGEBRA_COEFFICIENTS_H
