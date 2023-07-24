/* *************************************************************

Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurkó and Arend Janssen.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */

#ifndef DJC_COROPA_LIBALGEBRA_POLYLIEBASISH_SEEN
#define DJC_COROPA_LIBALGEBRA_POLYLIEBASISH_SEEN
#include "base_basis.h"
#include "basis.h"

namespace alg {

/// A basis for the polynomial Lie algebra, poly_lie.

/** A basis for the polynomial lie algebra.
The basis elements are vector fields of the of the form std::pair(direction,
monomial). The product is the Lie bracket of two vector fields.
*/

template<DEG n_letters, DEG max_degree>
class poly_lie_basis
{
public:
    /// The basis elements of poly_basis.
    typedef poly_basis POLYBASIS;
    /// The key of poly_basis (ie a monomial).
    typedef typename POLYBASIS::KEY POLYBASIS_KEY;
    /// A key is a pair of letter and monomial (ie a monomial in direction
    /// letter).
    typedef std::pair<LET, POLYBASIS_KEY> KEY;
    /// Polynomial algebra.

    /// The order in the MAP class reflects the degree
    struct KEY_LESS {
        bool inline operator()(const KEY& lhs, const KEY& rhs) const
        {
            return ((degree(lhs) < degree(rhs)) || ((degree(lhs) == degree(rhs)) && lhs < rhs));
        }
    };

public:
    // Property tags
    typedef alg::basis::with_degree<max_degree> degree_tag;
    typedef alg::basis::ordered<KEY_LESS> ordering_tag;

public:
    /// Default constructor. Empty basis.
    poly_lie_basis()
    {}

public:
    /// Turns a d/dx_i into a polynomial vector field by multiplying the empty
    /// monomial by d/dx_i.
    inline static KEY keyofletter(LET letter)
    {
        POLYBASIS empty;
        KEY result(letter, empty.empty_key);
        return result;
    }

    /// Returns the degree of the monomial in the pair (Let, monomial)
    inline static DEG degree(const KEY& k)
    {
        return POLYBASIS::degree(k.second);
    }

    /// Outputs a std::pair<poly_basis*, KEY> to an std::ostream.
    inline friend std::ostream& operator<<(std::ostream& os, const std::pair<poly_lie_basis*, KEY>& t)
    {
        POLYBASIS poly1;
        std::pair<POLYBASIS*, POLYBASIS_KEY> polypair;
        polypair.first = &poly1;
        polypair.second = t.second.second;
        os << "{" << polypair << "}"
           << "d/dx" << t.second.first << "}";
        return os;
    }
};

namespace vectors {

template<DEG n_letters, DEG max_degree, typename Field>
struct vector_type_selector<poly_lie_basis<n_letters, max_degree>, Field> {
    typedef poly_lie_basis<n_letters, max_degree> BASIS;
    typedef sparse_vector<BASIS, Field, std::map<typename BASIS::KEY, typename Field::S, typename BASIS::KEY_LESS>> type;
};

template<DEG n_letters, DEG max_degree, typename Field>
struct template_vector_type_selector<poly_lie_basis<n_letters, max_degree>, Field> {
    typedef poly_lie_basis<n_letters, max_degree> BASIS;

    template<typename B, typename C>
    using type = sparse_vector<B, C, std::map<typename B::KEY, typename C::S, typename B::KEY_LESS>>;
};

}// namespace vectors

namespace basis {

template<DEG Width, DEG Depth1, DEG Depth2>
struct related_to<poly_lie_basis<Width, Depth1>, poly_lie_basis<Width, Depth2>>
    : std::true_type {
};

}// namespace basis

}// namespace alg
#endif
