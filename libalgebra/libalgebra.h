/* *************************************************************

Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurkó and Arend Janssen.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */

//  libalgebra.h

// Include once wrapper
#ifndef DJC_COROPA_LIBALGEBRAH_SEEN
#define DJC_COROPA_LIBALGEBRAH_SEEN



#include <algorithm>
#include <cassert>
#include <cmath>
#include <deque>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <utility>
#include <vector>

#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
#include <boost/serialization/base_object.hpp>
#endif
#include <boost/thread/locks.hpp>
#include <boost/thread/recursive_mutex.hpp>
#include <boost/thread/shared_mutex.hpp>



/// The libalgebra namespace. A set of template classes for Algebras.
/**
   The SCA typename corresponds to a commutative ring with unit containing the
   rationals. Needed operators are * + - and explicit ctor from int type.

   The RAT typename corresponds to rationals built from SCA. Needed
   operators are + - * / and explicit ctor from SCA type and int type.
*/

#include "implementation_types.h"

#include "hall_set.h"
#include "detail/integer_maths.h"
#include "vectors.h"
#include "multiplication_helpers.h"

#include "algebra.h"

namespace alg {

/// Forward declaration of classes

/// Sparse vectors with default MAP typename from BASIS typename.
// template<class BASIS, class MAP = typename BASIS::MAP>
// class sparse_vector;
/// Generic Associative Algebra.
/*
template<typename Basis,
         typename Coeff,
         typename Multiplication,
         template<typename, typename, typename...> class VectorType = alg::vectors::template_vector_type_selector<Basis, Coeff>::template type,
         typename...>
class algebra;
*/

/// Generic Associative Algebra basis.
template<DEG n_letters, DEG max_degree = 0>
class tensor_basis;

/// Free Associative Algegra Basis. Concatenation product. Non commutative.
template<DEG n_letters, DEG max_degree = 0>
class free_tensor_basis;

/// Free Shuffle Associative Algebra Basis. Shuffle product. Commutative.
template<DEG n_letters, DEG max_degree = 0>
class shuffle_tensor_basis;


template <DEG Width, DEG Depth>
class traditional_free_tensor_multiplication;

template <DEG Width, DEG Depth, IDEG TileLetters=0, IDEG WriteCacheLetters=0>
class tiled_free_tensor_multiplication;

template <DEG Width, DEG Depth>
using free_tensor_multiplication = traditional_free_tensor_multiplication<Width, Depth>;

/// Free Associative Algebra.  Associative and non commutative.
template<typename Coeff, DEG n_letters, DEG max_degree = 0,
         template<typename, typename, typename...> class VectorType = vectors::template_vector_type_selector<free_tensor_basis<n_letters, max_degree>, Coeff>::template type,
         template <DEG, DEG> class FTMultiplication=free_tensor_multiplication,
         typename...>
class free_tensor;

/// Free Associative Shuffle Algebra.  Associative and Commutative.
template<typename Coeff, DEG n_letters, DEG max_degree = 0,
         template<typename, typename, typename...> class VectorType = vectors::template_vector_type_selector<shuffle_tensor_basis<n_letters, max_degree>, Coeff>::template type,
         typename...
         >
class shuffle_tensor;

/// Philip Hall Lie Basis.
template<DEG n_letters, DEG MaxDegree>
class hall_basis;

/// Free Lie Associative Algebra Basis.  Associative and non commutative.
template<DEG n_letters, DEG max_degree = 0>
class lie_basis;

/// Free Lie Associative Algebra.  Associative and non commutative.
template<typename Coeff, DEG n_letters, DEG max_degree = 0,
         template<typename, typename, typename...> class VectorType = vectors::template_vector_type_selector<lie_basis<n_letters, max_degree>, Coeff>::template type,
         typename...>
class lie;

/// Maps between Free Lie and Free Algebra elements.
template<typename Coeff, DEG n_letters, DEG max_degree = 0,
         typename Tensor = free_tensor<Coeff, n_letters, max_degree>,
         typename Lie = lie<Coeff, n_letters, max_degree>>
class maps;

/// Campbell-Baker-Hausdorff formulas.
template<typename Coeff, DEG n_letters, DEG max_degree,
         typename Tensor = free_tensor<Coeff, n_letters, max_degree>,
         typename Lie = lie<Coeff, n_letters, max_degree>>
class cbh;

/// Multivariate Polynomial Algebra Basis. Associative and Commutative.
class poly_basis;

/// Multivariate Polynomial Algebra.  Associative and Commutative.
template<typename Coeff,
         /*template <typename, typename, typename...> VectorType,*/
         typename...>
class poly;

/// II. Multivariate Polynomial Algebra pre Basis. Associative and Commutative
template<DEG n_letters, DEG max_degree = 0>
class monomial_basis;

/// II. Multivariate Polynomial Algebra Basis. Associative and Commutative
template<DEG n_letters, DEG max_degree = 0>
class free_monomial_basis;

/// II. Multivariate Polynomial Algebra   Associative and Commutative.
template<typename Coeff, DEG n_letters, DEG max_degree,
         template<typename, typename, typename...> class VectorType = vectors::template_vector_type_selector<free_monomial_basis<n_letters, max_degree>, Coeff>::template type,
         typename... Args>
class multi_polynomial;

/// III. Multivariate Polynomial Lie Algebra Basis. Associative and non
/// commutative
template<DEG n_letters, DEG max_degree = 0>
class poly_lie_basis;

/// III. Multivariate Polynomial Lie Algebra. Associative and non commutative
template<typename Coeff, DEG n_letters, DEG max_degree,
         template<typename, typename, typename...> class VectorType = vectors::template_vector_type_selector<poly_lie_basis<n_letters, max_degree>, Coeff>::template type,
         typename... Args>
class poly_lie;

}// namespace alg

#include "lie.h"
#include "lie_basis.h"
#include "monomial_basis.h"
#include "multi_polynomial.h"
#include "poly_basis.h"
#include "poly_lie.h"
#include "poly_lie_basis.h"
#include "polynomials.h"
#include "tensor.h"
#include "tensor_basis.h"
#include "utils.h"

#include "mpfloat_coefficients.h"
#include "rational_coefficients.h"
#include "functionals.h"
#include "multi_linear_operators.h"
#include "operators.h"
#include "alternative_multiplications.h"
#include "half_shuffle_tensor_basis.h"
#include "half_shuffle_tensor_multiplication.h"
#include "area_tensor_basis.h"
#include "area_tensor_multiplication.h"

// End of undeclaration of local macros.

// Include once wrapper
#endif// DJC_COROPA_LIBALGEBRAH_SEEN

// EOF.
