//
// Created by sam on 10/02/2021.
//

#ifndef LIBALGEBRAUNITTESTS_INTEGER_MATHS_H
#define LIBALGEBRAUNITTESTS_INTEGER_MATHS_H

namespace alg {
namespace integer_maths {

template<typename Base, typename Exp>
constexpr Base power(Base base, Exp exponent)
{
    return (exponent == 0) ? static_cast<Base>(1) : (exponent == 1) ? base
                                                                    : ((exponent % 2) ? base : static_cast<Base>(1)) * power(base, exponent / 2) * power(base, exponent / 2);
}

// Greatly simplify the computation of the mobius function by constexpr functions

template<typename Unsigned>
constexpr bool is_squarefree_impl(Unsigned N, Unsigned base)
{
    return (N % (base * base)) != 0 && (base * base >= N || is_squarefree_impl(N, base + 1));
}

template<typename Unsigned>
constexpr bool is_squarefree(Unsigned N)
{
    return is_squarefree_impl(N, static_cast<Unsigned>(2));
}

template<typename Unsigned>
constexpr typename std::make_signed<Unsigned>::type mobius(Unsigned);

// Mobius function is defined recursively
template<typename Unsigned>
constexpr typename std::make_signed<Unsigned>::type mobius_impl(Unsigned N, Unsigned divisor)
{
    using Int = typename std::make_signed<Unsigned>::type;
    return (divisor == N) ? static_cast<Int>(-1) : ((N % divisor == static_cast<Unsigned>(0)) ? mobius(divisor) * mobius(N / divisor)// Is a divisor, do product of mobius function on divisor and N/divisor
                                                                                              : mobius_impl(N, divisor + 1));        // not a divisor, increase divisor.
}

template<typename Unsigned>
constexpr typename std::make_signed<Unsigned>::type mobius(Unsigned N)
{
    using Int = typename std::make_signed<Unsigned>::type;
    return (N == static_cast<Unsigned>(1))
            ? static_cast<Int>(1)// mobius(1) = 1
            : (!is_squarefree(N) ? static_cast<Int>(0) : mobius_impl(N,
                                                                     static_cast<Unsigned>(2)));// mobius(as^2) = 0 for any a, s > 1;  otherwise recurse
}

template <typename ArgInt, typename BaseInt>
constexpr BaseInt logN(ArgInt arg, BaseInt base) noexcept
{
    return (arg < base) ? 0 : logN(arg / static_cast<ArgInt>(base), base) + 1;
}

template<typename Int, typename PowerInt>
constexpr Int sum_powers(Int base, PowerInt max_power) noexcept
{
    return (max_power == 0) ? Int(1) : (max_power == 1) ? Int(1) + base
                                                        : Int(1) + base * sum_powers(base, max_power - 1);
}

}// namespace integer_maths

template<typename Int>
constexpr Int mobius(Int n) noexcept;

namespace dtl {

template<typename Int>
constexpr bool is_squarefree_impl(Int n, Int base)
{
    return (n % (base * base)) != 0 && (base * base >= n || is_squarefree_impl(static_cast<Int>(n), static_cast<Int>(base + 1)));
}

template<typename Int>
constexpr Int mobius_impl(Int n, Int divisor) noexcept
{
    return (divisor == n) ? -1 : (n % divisor == 0) ? mobius(divisor) * mobius(n / divisor)
                                                    : mobius_impl(n, divisor + 1);
}

}// namespace dtl

template<typename Int>
constexpr bool is_squarefree(Int n) noexcept
{
    return dtl::is_squarefree_impl(n, static_cast<Int>(2));
}

template<typename Int>
constexpr Int mobius(Int n) noexcept
{
    return (n == 1) ? 1 : (!is_squarefree(n) ? 0 : dtl::mobius_impl(n, static_cast<Int>(2)));
}

template<typename Int, typename PowerInt>
constexpr Int power(Int base, PowerInt exponent) noexcept
{
    return (exponent == 0)    ? 1
            : (exponent == 1) ? base
                              : ((exponent % 2 != 0) ? base : 1) * power(base, exponent / 2) * power(base, exponent / 2);
}


template <typename Int, typename PowerInt>
constexpr Int sum_powers(Int base, PowerInt max_power) noexcept
{
    return (max_power == 0) ? 1 : (max_power == 1) ? 1 + base : 1 + base*sum_powers(base, max_power-1);
}



}// namespace alg

#endif// LIBALGEBRAUNITTESTS_INTEGER_MATHS_H
