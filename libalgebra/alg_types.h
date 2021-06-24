/* *************************************************************

Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai, 
Greg Gyurk� and Arend Janssen. 

Distributed under the terms of the GNU General Public License, 
Version 3. (See accompanying file License.txt)

************************************************************* */

#ifndef alg_types_h__
#define alg_types_h__

#include "libalgebra/libalgebra.h"
#include <addons/gmpwrapper.h>
#include <boost/multiprecision/gmp.hpp>
//#include "mtl/mtl.h"

#include <vector>

//#pragma warning(push)
//#pragma warning (disable : 800)
//#include "../addons/gmpwrapper.h"
//#pragma warning(pop)


enum coefficient_t
{
    Rational,
    DPReal,
    SPReal
};

enum vector_t
{
    Sparse,
    Dense,
    Hybrid
};




namespace {

template <vector_t VectorType>
struct vector_selector;

template <>
struct vector_selector<Sparse>
{
    template <typename Basis, typename Coeffs>
    struct selector
    {
        typedef LIBALGEBRA_DEFAULT_MAP_TYPE map_type;
        typedef alg::vectors::sparse_vector<Basis, Coeffs, map_type> type;
    };

};

template <>
struct vector_selector<Dense>
{
    template <typename Basis, typename Coeffs>
    struct selector
    {
        typedef std::vector<typename Coeffs::S> storage_type;
        typedef alg::vectors::dense_vector<Basis, Coeffs, storage_type> type;
    };

};

template <>
struct vector_selector<Hybrid>
{
    template <typename Basis, typename Coeffs>
    struct selector
    {
        typedef std::vector<typename Coeffs::S> storage_type;
        typedef LIBALGEBRA_DEFAULT_MAP_TYPE map_type;
        typedef alg::vectors::policy::basic_resize_policy policy_type;
        typedef alg::vectors::hybrid_vector<Basis, Coeffs, policy_type, storage_type, map_type> type;
    };

};


template<coefficient_t F>
struct Field;

template<>
struct Field<Rational>
{
    //typedef mpq_class S;
    //typedef mpq_class Q;
    typedef boost::multiprecision::number<boost::multiprecision::gmp_rational, boost::multiprecision::et_off> S;
    typedef S Q;
};

template<>
struct Field<DPReal>
{
    typedef double S;
    typedef double Q;
};

template<>
struct Field<SPReal>
{
    typedef float S;
    typedef float Q;
};

} // anon namespace

template<size_t D, size_t W, coefficient_t F = Rational, vector_t VectorType=Hybrid>
struct alg_types
{
    const static coefficient_t FIELD = F;
    typedef typename Field<F>::S S;
    typedef typename Field<F>::Q Q;
    typedef S SCA;
    typedef Q RAT;
    typedef alg::DEG DEG;
    typedef alg::LET LET;
    static const unsigned DEPTH = D;
    static const unsigned myDIM = W;
    static const unsigned ALPHABET_SIZE = W;
    typedef alg::poly<S, Q> MULTIPOLY1;
    typedef alg::free_tensor<S, Q, ALPHABET_SIZE, DEPTH> TENSOR;
    typedef alg::lie<S, Q, ALPHABET_SIZE, DEPTH> LIE;
    typedef alg::maps<S, Q, ALPHABET_SIZE, DEPTH> MAPS;
    typedef alg::cbh<S, Q, ALPHABET_SIZE, DEPTH> CBH;
    typedef alg::poly_lie<S, Q, ALPHABET_SIZE, DEPTH> POLYLIE;
    typedef alg::multi_polynomial<S, Q, ALPHABET_SIZE, DEPTH> MULTIPOLY;
    //typedef mtl::dense1D<RAT> mtlVector;
    //typedef typename mtl::matrix<RAT, mtl::rectangle<>, mtl::dense<>, mtl::row_major>::type mtlMatrix;
    //typedef typename mtl::matrix<RAT, mtl::diagonal<>, mtl::packed<>, mtl::row_major>::type mtlDiagMat;
};



////  alg_types.h : provides an interface to and sets consistent sets of basic algebraic types
//
//
//#ifndef alg_types_h__
//#define alg_types_h__
//
//#include "libalgebra.h"
//
//#pragma warning(push)
//#pragma warning (disable : 800)
//#include "../addons/gmpwrapper.h"
//#pragma warning(pop)
//
//
//enum coefficient_t
//{
//	Rational,
//	DPReal,
//	SPReal
//};
//
//namespace
//{
//
//template <coefficient_t F>
//struct Field;
//
//template<>
//struct Field<Rational>
//{
//	typedef mpq_class S;
//	typedef mpq_class Q;
//};
//
//template<>
//struct Field<DPReal>
//{
//	typedef double S;
//	typedef double Q;
//};
//
//template<>
//struct Field<SPReal>
//{
//	typedef float S;
//	typedef float Q;
//};
//
//} // anon namespace
//
//template <unsigned D, unsigned W, coefficient_t F = Rational>
//struct alg_types : Field<F>
//{
//	typedef typename Field<F>::S S;
//	typedef typename Field<F>::Q Q;
//	typedef S SCA;
//	typedef Q RAT;
//	static const unsigned DEPTH = D;
//	static const unsigned myDIM = W;
//	static const unsigned ALPHABET_SIZE = W;
//	typedef alg::poly<S, Q> MULTIPOLY1;
//	typedef alg::free_tensor<S, Q, ALPHABET_SIZE, DEPTH> TENSOR;
//	typedef alg::lie<S, Q, ALPHABET_SIZE, DEPTH> LIE;
//	typedef alg::maps<S, Q, ALPHABET_SIZE, DEPTH> MAPS;
//	typedef alg::cbh<S, Q, ALPHABET_SIZE, DEPTH> CBH;
//	typedef alg::multi_polynomial<S, Q, ALPHABET_SIZE, DEPTH> MULTIPOLY;
//};

#endif // alg_types_h__
