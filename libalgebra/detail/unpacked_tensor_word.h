//
// Created by user on 05/12/22.
//

#ifndef LIBALGEBRA_LIBALGEBRA_DETAIL_UNPACKED_TENSOR_WORD_H_
#define LIBALGEBRA_LIBALGEBRA_DETAIL_UNPACKED_TENSOR_WORD_H_

#include "smallest_int_type.h"
#include <cstdint>
#include <cassert>
#include <algorithm>

namespace alg {

template <unsigned Width, unsigned Depth>
class unpacked_tensor_word
{
    using letter_type = typename dtl::smallest_containing_int<Width>::type;
    using index_type = std::ptrdiff_t;
    using size_type = std::size_t;

    letter_type m_bits[Depth] = {};
    unsigned m_degree;

public:

    template <typename Int>
    void reset(Int index) noexcept
    {
        m_degree = unsigned(index);
        std::fill(m_bits, m_bits+Depth, letter_type(0));
    }

    constexpr unpacked_tensor_word& operator++() noexcept
    {
        letter_type carry = 1;
        int pos = m_degree;
        do {
            --pos;
            m_bits[pos] += carry;
            if (m_bits[pos] >= Width) {
                m_bits[pos] = letter_type(0);
            } else {
                carry = 0;
            }
        } while (carry > 0 && pos > 0);
        m_degree += carry;
        return *this;
    }


    template <typename IndexType=index_type>
    constexpr IndexType to_index() const noexcept
    {
        IndexType result = 0;
        for (auto i=1; i<=m_degree; ++i) {
            result *= IndexType(Width);
            result += m_bits[m_degree-i];
        }
        return result;
    }

    template <typename IndexType=index_type>
    constexpr IndexType to_reverse_index() const noexcept
    {
        IndexType result = 0;
        for (auto i=0; i<m_degree; ++i) {
            result *= IndexType(Width);
            result += m_bits[i];
        }
        return result;


    }

    template <typename IndexType=index_type, typename Int>
    constexpr IndexType split_left_index(Int left_letters) const noexcept
    {
        IndexType result = 0;
        for (auto i=0; i<left_letters; ++i) {
            result *= IndexType(Width);
            result += m_bits[i];
        }
        return result;
    }

    template <typename IndexType=index_type, typename Int>
    constexpr IndexType split_left_reverse_index(Int left_letters) const noexcept
    {
        IndexType result = 0;
        for (auto i=1; i<=left_letters; ++i) {
            result *= Width;
            result += m_bits[left_letters-i];
        }
        return result;
    }

    template <typename IndexType=index_type, typename Int>
    constexpr IndexType split_right_index(Int left_letters) const noexcept
    {
        IndexType result = 0;
        for (auto i=left_letters; i<m_degree; ++i) {
            result *= Width;
            result += m_bits[i];
        }
        return result;
    }


};


} // namespace alg

#endif//LIBALGEBRA_LIBALGEBRA_DETAIL_UNPACKED_TENSOR_WORD_H_
