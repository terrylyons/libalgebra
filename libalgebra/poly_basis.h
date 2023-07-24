/* *************************************************************

Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurkó and Arend Janssen.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */

//  poly_basis.h

// Include once wrapper
#ifndef DJC_COROPA_LIBALGEBRA_POLYBASISH_SEEN
#define DJC_COROPA_LIBALGEBRA_POLYBASISH_SEEN

#include "basis.h"
#include "implementation_types.h"
#include "vectors.h"

#include <iosfwd>
#include <map>

namespace alg {

/// A polynomial basis class.
/**
This is the basis used to implement the polynomial class as a
specialisation of the algebra class. A key implements a monomial in several
variables, with scalar coefficient one. Each variable corresponds to a
letter. The product of monomial, provided by the member function prod(), is
the standard commutative product of monomials. The type LET is used to
number the variables (i.e. letters). The type SCA is used to evaluate the
monomial. In the current implementation, no key is stored in memory, only
the member functions are provided.

An empty monomial corresponds to a constant term, i.e. a scalar SCA.
*/

class poly_basis
{
public:
    typedef alg::LET LET;
    /// A key is a map from letters to degrees (i.e. a monomial of letters).
    typedef std::map<LET, DEG> KEY;
    /// A default key corresponds to the monomial of degree 0.
    const KEY empty_key;
    /// The rationals.

    struct KEY_LESS;
    /// The MAP type.

public:
    // Property tags
    typedef alg::basis::without_degree degree_tag;
    typedef alg::basis::ordered<KEY_LESS> ordering_tag;

public:
    /// Default constructor. Empty basis.
    poly_basis(void)
    {}

public:
    // tjl the code below is strange -
    // the result of such an evaluation should be a scalar times a key of
    // unevaluated variables
    /// Evaluate the key from a vector of scalars.
    template<typename Coeff>
    typename Coeff::SCA eval_key(const KEY& k, const std::map<LET, typename Coeff::SCA>& values) const
    {
        typename Coeff::SCA result(Coeff::one);
        typename KEY::const_iterator it;
        for (it = k.begin(); it != k.end(); ++it) {
            if (it->second > 0) {
                typename std::map<LET, typename Coeff::SCA>::const_iterator iit;
                iit = values.find(it->first);
                if (iit != values.end()) {
                    for (DEG j = 1; j <= it->second; ++j) {
                        Coeff::mul_inplace(result, iit->second);
                    }
                }
                else {
                    throw std::invalid_argument("not all variables have values!");
                }
            }
        }

        return result;
    }

    /// Returns the key (monomial) corresponding to a letter (variable).
    inline KEY keyofletter(LET letter) const
    {
        KEY result;
        result[letter] = +1;
        return result;
    }

    /// Returns the degree of a monomial
    inline static DEG degree(const KEY& k)
    {
        DEG result(0);
        typename KEY::const_iterator it;
        for (it = k.begin(); it != k.end(); ++it) {
            result += DEG(it->second);
        }
        return result;
    }

    struct KEY_LESS {
        bool inline operator()(const KEY& lhs, const KEY& rhs) const
        {
            return ((degree(lhs) < degree(rhs)) || ((degree(lhs) == degree(rhs)) && lhs < rhs));
        }
    };

    /// Returns the value of the smallest key in the basis.
    inline KEY begin(void) const
    {
        return empty_key;
    }
    ///// Returns the key next the biggest key of the basis.
    //// this doesn't make sense without a maximum degree.
    // inline KEY end(void) const
    //{
    // KEY result; // empty key.
    // result.push_back(0); // invalid key.
    // return result;
    // }
    ///// Returns the key next a given key in the basis.
    //// We need an ordering for this to work
    // inline KEY nextkey(const KEY& k) const
    //{
    // KEY::size_type i;
    // for (i = k.size()-1; i >= 0; --i)
    // if (k[i]<n_letters) { KEY result(k); result[i] += 1; return result; }
    // return end();
    // }

    /// Outputs a std::pair<poly_basis*, KEY> to an std::ostream.
    inline friend std::ostream& operator<<(std::ostream& os, const std::pair<poly_basis*, KEY>& t)
    {
        bool first(true);
        typename KEY::const_iterator it;
        for (it = t.second.begin(); it != t.second.end(); ++it) {
            if (it->second > 0) {
                if (!first) {
                    os << ' ';
                }
                else {
                    first = false;
                }
                os << "x" << it->first;
                if (it->second > 1) {
                    os << '^' << it->second;
                }
            }
        }
        return os;
    }
};

namespace vectors {

template<typename Coeff>
struct vector_type_selector<poly_basis, Coeff> {
    typedef poly_basis BASIS;
    typedef sparse_vector<BASIS, Coeff, std::map<typename BASIS::KEY, typename Coeff::S, typename BASIS::KEY_LESS>> type;
};

template<typename Coeff>
struct template_vector_type_selector<poly_basis, Coeff> {
    typedef poly_basis BASIS;
    template<typename B, typename C>
    using type = sparse_vector<B, C, std::map<typename B::KEY, typename C::S, typename B::KEY_LESS>>;
};

}// namespace vectors

}// namespace alg

// Include once wrapper
#endif// DJC_COROPA_LIBALGEBRA_POLYBASISH_SEEN

// EOF.
