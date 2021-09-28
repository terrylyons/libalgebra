//
// Created by sam on 10/02/2021.
//

#ifndef LIBALGEBRAUNITTESTS_INTEGER_MATHS_H
#define LIBALGEBRAUNITTESTS_INTEGER_MATHS_H

namespace alg {
namespace integer_maths {

template <typename Base, typename Exp>
constexpr Base power(Base base, Exp exponent)
{
    return (exponent==0) ? static_cast<Base>(1) : (exponent==1) ? base :
                               ((exponent%2) ? base : static_cast<Base>(1))*power(base, exponent/2)*power(base, exponent/2);
}


// Greatly simplify the computation of the mobius function by constexpr functions

template <typename Unsigned>
constexpr bool is_squarefree_impl(Unsigned N, Unsigned base)
{
    return (N%(base*base))!=0 && (base*base>=N || is_squarefree_impl(N, base+1));
}

template <typename Unsigned>
constexpr bool is_squarefree(Unsigned N)
{
    return is_squarefree_impl(N, static_cast<Unsigned>(2));
}

template <typename Unsigned>
constexpr std::make_signed_t<Unsigned> mobius(Unsigned);

// Mobius function is defined recursively
template <typename Unsigned>
constexpr std::make_signed_t<Unsigned> mobius_impl(Unsigned N, Unsigned divisor)
{
    using Int = std::make_signed_t<Unsigned>;
    return (divisor==N) ? static_cast<Int>(-1) : ((N%divisor==static_cast<Unsigned>(0))
        ? mobius(divisor)*mobius(N/divisor) // Is a divisor, do product of mobius function on divisor and N/divisor
        : mobius_impl(N, divisor+1)); // not a divisor, increase divisor.
}

template <typename Unsigned>
constexpr std::make_signed_t<Unsigned> mobius(Unsigned N)
{
    using Int = std::make_signed_t<Unsigned>;
    return (N==static_cast<Unsigned>(1))
        ? static_cast<Int>(1) // mobius(1) = 1
        : (!is_squarefree(N) ? static_cast<Int>(0) : mobius_impl(N, static_cast<Unsigned>(2))); // mobius(as^2) = 0 for any a, s > 1;  otherwise recurse
}





} // namespace integer_maths

template <typename Int>
constexpr Int mobius(Int n) noexcept;



namespace dtl {

    template<typename Int>
    constexpr bool is_squarefree_impl(Int n, Int base)
    {
        return (n%(base*base))!=0 &&
            (base*base>=n || is_squarefree_impl(static_cast<Int>(n), static_cast<Int>(base+1)));
    }

    template <typename Int>
    constexpr Int mobius_impl(Int n, Int divisor) noexcept
    {
        return (divisor == n) ? -1 :
               (n % divisor == 0) ? mobius(divisor) * mobius(n / divisor)
                                  : mobius_impl(n, divisor+1);
    }

}

template<typename Int>
constexpr bool is_squarefree(Int n) noexcept
{
    return dtl::is_squarefree_impl(n, static_cast<Int>(2));
}

template <typename Int>
constexpr Int mobius(Int n) noexcept
{
    return (n == 1) ? 1 : (!is_squarefree(n) ? 0 : dtl::mobius_impl(n, static_cast<Int>(2)));
}


template <typename Int, typename PowerInt>
constexpr Int power(Int base, PowerInt exponent) noexcept
{
    return (exponent == 0) ? 1
            : (exponent == 1) ? base :
              ((exponent % 2 != 0) ? base : 1) * power(base, exponent / 2) * power(base, exponent /2);
}


} // namespace alg

#endif // LIBALGEBRAUNITTESTS_INTEGER_MATHS_H
