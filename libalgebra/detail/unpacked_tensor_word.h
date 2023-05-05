//
// Created by user on 05/12/22.
//

#ifndef LIBALGEBRA_LIBALGEBRA_DETAIL_UNPACKED_TENSOR_WORD_H_
#define LIBALGEBRA_LIBALGEBRA_DETAIL_UNPACKED_TENSOR_WORD_H_

#include "smallest_int_type.h"
#include <cstdint>
#include <cassert>
#include <array>
#include <algorithm>

namespace alg {

template <unsigned Width, unsigned Depth>
class unpacked_tensor_word
{
    using letter_type = typename dtl::smallest_containing_int<Width>::type;
    using index_type = std::ptrdiff_t;
    using size_type = std::size_t;

    std::array<letter_type, Depth> m_bits {};
    unsigned m_degree;

public:

    constexpr unsigned degree() const noexcept { return m_degree; }

    template <typename Int>
    constexpr letter_type operator[](Int idx) const noexcept
    { return m_bits[idx]; }


    template <typename Int>
    void reset(Int index) noexcept
    {
        m_degree = unsigned(index);
        std::fill(m_bits.begin(), m_bits.end(), letter_type(0));
    }

    constexpr unpacked_tensor_word& operator++() noexcept
    {
        letter_type carry = 1;
        if (m_degree > 0) {
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
        }
        m_degree += carry;
        return *this;
    }


    template <typename IndexType=index_type>
    constexpr IndexType to_index() const noexcept
    {
        IndexType result = 0;
        for (unsigned i=0; i<m_degree; ++i) {
            result *= IndexType(Width);
            result += m_bits[i];
        }
        return result;
    }

    template <typename IndexType=index_type>
    constexpr IndexType to_reverse_index() const noexcept
    {
        IndexType result = 0;
        for (unsigned i=1; i<=m_degree; ++i) {
            result *= IndexType(Width);
            result += m_bits[m_degree-i];
        }
        return result;


    }

    template <typename IndexType=index_type, typename Int>
    constexpr IndexType split_left_index(Int left_letters) const noexcept
    {
        IndexType result = 0;
        for (unsigned i=0; i<left_letters; ++i) {
            result *= IndexType(Width);
            result += m_bits[i];
        }
        return result;
    }

    template <typename IndexType=index_type, typename Int>
    constexpr IndexType split_left_reverse_index(Int left_letters) const noexcept
    {
        IndexType result = 0;
        for (Int i=1; i<=left_letters; ++i) {
            result *= IndexType(Width);
            result += m_bits[left_letters-i];
        }
        return result;
    }

    template <typename IndexType=index_type, typename Int>
    constexpr IndexType split_right_index(Int left_letters) const noexcept
    {
        IndexType result = 0;
        for (unsigned i=left_letters; i<m_degree; ++i) {
            result *= IndexType(Width);
            result += m_bits[i];
        }
        return result;
    }

    template <typename IndexType=index_type, typename Int>
    constexpr IndexType split_right_reverse_index(Int left_letters) const noexcept
    {
        IndexType result = 0;
        assert(left_letters >= 0 && left_letters <= m_degree);
        for (unsigned i=m_degree; i > left_letters; --i) {
            result *= IndexType(Width);
            result += IndexType(m_bits[i-1]);
        }
        return result;
    }



};

template <unsigned Depth>
struct unpacked_tensor_word<1U, Depth> {
    using letter_type = unsigned char;
    using index_type = std::ptrdiff_t;
    using size_type = std::size_t;

    unsigned m_degree;

public:
    constexpr unsigned degree() const noexcept { return m_degree; }

    template<typename Int>
    constexpr letter_type operator[](Int idx) const noexcept
    {
        return 0;
    }

    template<typename Int>
    void reset(Int index) noexcept{
    }

    constexpr unpacked_tensor_word& operator++() noexcept
    {
        m_degree += 1;
        return *this;
    }

    template<typename IndexType = index_type>
    constexpr IndexType to_index() const noexcept
    {
        return m_degree;
    }

    template<typename IndexType = index_type>
    constexpr IndexType to_reverse_index() const noexcept
    {
        return m_degree;
    }

    template<typename IndexType = index_type, typename Int>
    constexpr IndexType split_left_index(Int left_letters) const noexcept
    {
        return std::min(m_degree, static_cast<unsigned>(left_letters));
    }

    template<typename IndexType = index_type, typename Int>
    constexpr IndexType split_left_reverse_index(Int left_letters) const noexcept
    {
        return std::min(m_degree, static_cast<unsigned>(left_letters));
    }

    template<typename IndexType = index_type, typename Int>
    constexpr IndexType split_right_index(Int left_letters) const noexcept
    {
        return left_letters > m_degree ? 0 : m_degree - left_letters;
    }

    template<typename IndexType = index_type, typename Int>
    constexpr IndexType split_right_reverse_index(Int left_letters) const noexcept
    {
        return left_letters > m_degree ? 0 : m_degree - left_letters;
    }

};

template <unsigned Width>
struct unpacked_tensor_word<Width, 0U> {
    using letter_type = unsigned char;
    using index_type = std::ptrdiff_t;
    using size_type = std::size_t;

public:

    constexpr unsigned degree() const noexcept { return 0; }

    template<typename Int>
    constexpr letter_type operator[](Int idx) const noexcept
    {
        return 0;
    }

    template<typename Int>
    void reset(Int index) noexcept
    {
    }

    constexpr unpacked_tensor_word& operator++() noexcept
    {
        return *this;
    }

    template<typename IndexType = index_type>
    constexpr IndexType to_index() const noexcept
    {
        return 0;
    }

    template<typename IndexType = index_type>
    constexpr IndexType to_reverse_index() const noexcept
    {
        return 0;
    }

    template<typename IndexType = index_type, typename Int>
    constexpr IndexType split_left_index(Int left_letters) const noexcept
    {
        return 0;
    }

    template<typename IndexType = index_type, typename Int>
    constexpr IndexType split_left_reverse_index(Int left_letters) const noexcept
    {
        return 0;
    }

    template<typename IndexType = index_type, typename Int>
    constexpr IndexType split_right_index(Int left_letters) const noexcept
    {
        return 0;
    }

    template<typename IndexType = index_type, typename Int>
    constexpr IndexType split_right_reverse_index(Int left_letters) const noexcept
    {
        return 0;
    }
};


} // namespace alg

#endif//LIBALGEBRA_LIBALGEBRA_DETAIL_UNPACKED_TENSOR_WORD_H_
