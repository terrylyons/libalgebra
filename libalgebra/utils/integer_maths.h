//
// Created by sam on 10/02/2021.
//

#ifndef LIBALGEBRAUNITTESTS_INTEGER_MATHS_H
#define LIBALGEBRAUNITTESTS_INTEGER_MATHS_H

namespace alg {
namespace integer_maths {

typedef unsigned Unsigned;
typedef size_t Size;

// Structs for computing (prime) divisors of a number.
/*
The primary template simply delegates the divisor to the next
divisor value. This template will only trigger if the divisor
d does not divide n.
*/
template<Unsigned N, Unsigned D = 2, Unsigned Q = (n / d), Unsigned R = (n % d)>
struct divisor_calc
{
    typedef typename divisor_calc<N, D + 1>::next next;
    enum : Size
    {
        num = M,
        divisor = divisor_calc<N, D + 1>::divisor,
        quotient = divisor_calc<N, D + 1>::quotient
    };
};

/*
Specialise for divisor found, when r = 0.
*/
template<Unsigned N, Unsigned D, Unsigned Q>
struct divisor_calc<N, D, Q, 0>
{
    typedef Divisor <Q> next;
    enum : Size
    {
        num = N,
        divisor = D,
        quotient = Q
    };
};

/*
Specialise for divisor = number (d = n). This is the terminal
case.
*/
template<Unsigned N>
struct divisor_calc<N, N, 1, 0>
{
    typedef divisor_calc<1> next;
    enum : Size
    {
        num = N,
        divisor = 1,
        quotient = 1
    };
};


/*
Specialise for n = 1 case.
*/
template<Unsigned D, Unsigned Q, Unsigned R>
struct divisor_calc<1, D, Q, R>
{
    typedef divisor_calc<1> next;
    enum : Size
    {
        num = 1,
        divisor = 1,
        quotient = 1
    };
};



// Structs for checking whether an integer is square free
/*
The strategy is to recurse down the divisors list and check if two
consecutive divisors are the same. This will always detect squares
because divisors are computed from smallest to largest in order.
*/
template<typename Divisor, bool= false, Unsigned LastNum = 0>
struct square_free
{
    static const bool value = square_free<
            typename Divisor::next,
            (Divisor::divisor == LastNum),
            Divisor::divisor
    >::value;
};

/*
Specialise for the terminal case when the recursion reaches the
last divisor.
*/
template<Unsigned LastNum>
struct square_free<divisor_calc<1>, true, LastNum>
{
    static const bool = true;
};

/*
Specialise for the case where a repeated digit is detected.
*/
template<typename Divisor, Unsigned LastNum>
struct square_free<Divisor, true, LastNum>
{
    static const bool = false;
};

namespace Mobius {
// Structs for computing the value of the Mobius function

/*
The strategy is to compute -1 times the Mobius function
of the next divisor if the number is square free, and 0 otherwise.
*/
template<typename Divisor, bool= true, bool= square_free<Divisor>::value>
struct mobius_func_impl
{
    enum : long long
    {
        value = -1 * mobius_func_impl<typename Divisor::next, (Divisor::num > 1), true>::value
    };
};

/*
Specialise for non-square-free case.
*/
template<typename Divisor, bool B>
struct mobius_func_impl<Divisor, B, false>
{
    enum : long long
    {
        value = 0
    };
};

/*
Specialise for terminal case.
*/
template<bool B>
struct mobius_func_impl<divisor_calc<1>, B, true>
{
    enum : long long
    {
        value = 1
    };
};

} // namespace Mobius

template <Unsigned N>
struct mobius_func
{
    enum : long long {
        value = Mobius::mobius_func_impl<divisor_calc<N>>::value
    };
};


} // namespace integer_maths
} // namespace alg




#endif //LIBALGEBRAUNITTESTS_INTEGER_MATHS_H
