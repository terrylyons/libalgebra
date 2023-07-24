/* *************************************************************

Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurkó and Arend Janssen.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */

#pragma once
#ifndef ConstPower_h__
#define ConstPower_h__

namespace alg {
/// A template for constructing integer constants

/// ConstPower< arg, exp>::ans is the constant integer value of arg^exp
template<unsigned long long arg, DEG exp>
struct ConstPower {
    enum : long long
    {
        enum1_force_long_long = 9223372036854775807ULL,
        ans = ConstPower<arg, exp / 2>::ans * ConstPower<arg, exp / 2>::ans * ((exp % 2 == 0) ? 1ULL : arg),

    };
};

template<unsigned long long arg>
struct ConstPower<arg, 0> {
    enum : long long
    {
        enum1_force_long_long = 9223372036854775807ULL,
        ans = 1ULL,
    };
};

/// Test of ConstPower Template

/// TestConstPower<3,5> tests the ConstPower Template for some the powers x^y
/// with x < 4 and y < 16
template<unsigned long long arg, DEG exp>
struct TestConstPower {
    enum
    {
        intermediate = ConstPower<arg, exp>::ans * TestConstPower<arg, exp - 1>::intermediate,
        ans = (intermediate == ConstPower<arg, (ConstPower<exp, 2>::ans + exp) / 2>::ans) && TestConstPower<arg - 1, exp>::ans,
        enum1_force_long_long = 9223372036854775807ULL
    };
};

template<unsigned long long arg>
struct TestConstPower<arg, 0> {
    enum
    {
        intermediate = ConstPower<arg, 0>::ans,
        ans = (intermediate == 1ULL),
        enum1_force_long_long = 9223372036854775807ULL
    };
};

template<DEG exp>
struct TestConstPower<1, exp> {
    enum
    {
        intermediate = ConstPower<1, exp>::ans,
        ans = (intermediate == 1)
    };
};

template<>
struct TestConstPower<1, 0> {
    enum
    {
        intermediate = ConstPower<1, 0>::ans,
        ans = (intermediate == 1)
    };
};

}// namespace alg
#endif// ConstPower_h__
