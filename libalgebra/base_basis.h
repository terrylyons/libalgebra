/* *************************************************************

Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurkó and Arend Janssen.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */

#pragma once
#ifndef LIBALGEBRA_BASE_BASIS_H__
#define LIBALGEBRA_BASE_BASIS_H__

namespace alg {

/// Basis traits

enum basis_product_trait
{
    Without_Degree,
    With_Degree,
    No_Product
};

/// The basis' properties

/// This structure stores the number of letters: NO_LETTERS, maximum degree:
/// MAX_DEGREE and the product_trait, either Without_Degree, With_Degree or
/// No_Product.
template<basis_product_trait trait = Without_Degree, DEG no_letters = 0, DEG max_degree = 0>
struct base_basis {

    /// The number of letters used to generate the algebra
    /// zero if there is no well defined finite generating set
    static constexpr DEG NO_LETTERS = no_letters;
    /// The trait of the product; either Without_Degree, With_Degree or No_Product
    static constexpr basis_product_trait PRODUCT_TYPE = trait;
    /// The maximum degree
    static constexpr DEG MAX_DEGREE = (trait == With_Degree) ? max_degree : DEG(0);
};

}// namespace alg

#endif// LIBALGEBRA_BASE_BASIS_H__
