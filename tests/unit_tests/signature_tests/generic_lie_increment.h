//
// Created by sam on 11/03/2021.
//

#ifndef LIBALGEBRAUNITTESTS_GENERIC_LIE_INCREMENT_H
#define LIBALGEBRAUNITTESTS_GENERIC_LIE_INCREMENT_H

#include <vector>
#include <boost/array.hpp>

#include "generic_coefficient.h"

template <unsigned Width, typename Integer>
class generic_lie_increment
{
    boost::array<generic_coefficient<Integer>, Width> m_data;
public:

    generic_lie_increment() : m_data(0)
    {
    }

    generic_lie_increment(const generic_lie_increment& other)
        : m_data(other.m_data)
    {}

    explicit
    generic_lie_increment(std::vector<generic_coefficient<Integer> > incr)
        : m_data()
    {
        assert(incr.size() == Width);
        for (size_t i=0; i<Width; ++i) {
            m_data[i] = incr[i];
        }
    }

    std::size_t size() const
    {
        return m_data.size();
    }

    template <typename Lie>
    Lie to_lie() const
    {
        typedef typename Lie::KEY Key; // size_t
        typedef typename Lie::SCALAR Scalar;
        Lie result;
        for (Key i=0; i<size(); ++i) {
            // zero is not a valid Lie key.
            result.add_scal_prod(i+1, Scalar(m_data[i]));
        }
        return result;
    }

    generic_coefficient<Integer>& operator[](std::size_t idx)
    {
        assert(idx < size());
        return m_data[idx];
    }

};




#endif //LIBALGEBRAUNITTESTS_GENERIC_LIE_INCREMENT_H
