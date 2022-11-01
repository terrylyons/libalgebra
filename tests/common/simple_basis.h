//
// Created by sam on 02/02/2021.
//

#ifndef LIBALGEBRAUNITTESTS_SIMPLE_BASIS_H
#define LIBALGEBRAUNITTESTS_SIMPLE_BASIS_H
#include <libalgebra/libalgebra.h>
#include <libalgebra/basis.h>

#include <iostream>
#include <utility>
#include <functional>

using alg::DIMN;
using alg::DEG;

template <unsigned D, typename R>
class SimpleIntegerBasis
{
public:
    static const unsigned dimension = D;
    typedef unsigned KEY;
    typedef R RATIONAL;

public:
    // Property tags
    typedef alg::basis::without_degree degree_tag;
    typedef alg::basis::ordered<std::less<KEY> > ordering_tag;


    friend std::ostream& operator<<(std::ostream& os,
            const std::pair<SimpleIntegerBasis*, KEY> arg)
    {
        return (os << arg.second);
    }

    DIMN key_to_index(const KEY k) const
    {
        return static_cast<DIMN>(k);
    }

    KEY index_to_key(const DIMN idx) const
    {
        return static_cast<KEY>(idx);
    }

    KEY nextkey(const KEY& k) const
    { return k + 1; }

    DIMN max_dimension()
    { return dimension; }

};




#endif //LIBALGEBRAUNITTESTS_SIMPLE_BASIS_H
