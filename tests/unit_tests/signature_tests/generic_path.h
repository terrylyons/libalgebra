//
// Created by sam on 11/03/2021.
//

#ifndef LIBALGEBRAUNITTESTS_GENERIC_PATH_H
#define LIBALGEBRAUNITTESTS_GENERIC_PATH_H

#include <vector>
#include <utility>
#include <functional>
#include <stdint.h>

#include "../../common/memfile.h"

#include "generic_lie_increment.h"

template <unsigned Width, typename Integer=int32_t>
class generic_path
{
    std::vector<generic_lie_increment<Width, Integer> > m_increments;
public:

    generic_path() : m_increments()
    {}

    generic_path(const generic_path& other) : m_increments(other.m_increments)
    {}

    generic_path(Integer* vals, size_t length) : m_increments()
    {
        m_increments.reserve(length);
        std::vector<generic_coefficient<Integer> > tmp;
        tmp.reserve(Width);

        for (int i=0; i < length; ++i) {
            for (int j=0; j<Width; ++j) {
                tmp.push_back(generic_coefficient<Integer>(vals[i*2*Width + 2*j], vals[i*2*Width + 2*j + 1]));
            }
            m_increments.push_back(generic_lie_increment<Width, Integer>(tmp));
            tmp.clear();
        }

    }

    explicit
    generic_path(std::vector<generic_lie_increment<Width, Integer> > increments)
        : m_increments(increments)
    {
    }

    size_t length() const
    {
        return m_increments.size();
    }

    template <typename Framework>
    typename Framework::TENSOR signature(size_t start_increment=0, size_t end_increment=-1) const
    {
        typedef typename Framework::TENSOR Tensor;
        typedef typename Framework::LIE Lie;
        typename Framework::MAPS maps;

        size_t end = std::min(end_increment, m_increments.size());
        assert(start_increment >= 0);
        assert(end <= m_increments.size());

        Tensor result(typename Tensor::SCALAR(1));
        for (size_t i=start_increment; i<end; ++i) {
            result *= exp(maps.l2t(m_increments[i].template to_lie<Lie>()));
            //result.fmexp_inplace(maps.l2t(m_increments[i].template to_lie<Lie>()));
        }
        return result;
    }

    template <typename Framework>
    typename Framework::LIE log_signature(size_t start_increment=0, size_t end_increment=-1) const
    {
        typename Framework::MAPS maps;
        return maps.t2l(log(signature<Framework>(start_increment, end_increment)));
    }

};





#endif //LIBALGEBRAUNITTESTS_GENERIC_PATH_H
