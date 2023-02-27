//
// Created by user on 18/10/22.
//

#ifndef LIBALGEBRA_LIBALGEBRA_DETAIL_REVERSING_PERMUTATION_H_
#define LIBALGEBRA_LIBALGEBRA_DETAIL_REVERSING_PERMUTATION_H_

#include <cstdint>

namespace alg { namespace dtl {

template<unsigned Width, unsigned Level>
struct reversing_permutation {
    using size_type = size_t;

    static constexpr size_type factor = power(Width, Level - 1);


    static constexpr size_type first_letter(size_type idx)
    {
        return idx / factor;
    }

    static constexpr size_type last_letter(size_type idx)
    {
        return idx % Width;
    }

    static constexpr size_type middle_word(size_type idx)
    {
        /*
             * Writing idx = l_2*Width^{Level-1} + index(middle_word)*Width + l1
             * we can rearrange to get index(middle_word) = (idx - l1 - l2*Width^{Level-1})/Width.
             * Although since l1 < Width, we can ignore it since floor division will take care of it.
             * The brackets on the right hand side can then be realised as idx % Width^{Level-1}
             */
        return (idx % factor) / Width;
    }

    static constexpr size_type permute_idx(size_type idx)
    {
        static_assert(Level - 2 > 0, "Level must be at least 3 in this specialisation");
        using next = reversing_permutation<Width, Level - 2>;

        constexpr size_type shift = power(Width, Level - 1);
        return last_letter(idx) * shift + next::permute_idx(middle_word(idx)) * Width + first_letter(idx);
    }

    template<typename T>
    void operator()(T* __restrict tile) const noexcept
    {
        for (size_type i = 0; i < power(Width, Level); ++i) {
            auto j = permute_idx(i);
            if (j > i) {
                std::swap(tile[i], tile[j]);
            }
        }
    }

};

template<unsigned Width>
struct reversing_permutation<Width, 2> {
    using size_type = size_t;
    static const unsigned Level = 2;

    static constexpr size_type first_letter(size_type idx)
    {
        return idx / Width;
    }

    static constexpr size_type last_letter(size_type idx)
    {
        return idx % Width;
    }

    static constexpr size_type permute_idx(size_type idx)
    {
        return last_letter(idx) * Width + first_letter(idx);
    }

    /// Operate inplace on a single tile
    template<typename T>
    void operator()(T* __restrict tile) const noexcept
    {
        for (size_type i = 0; i < power(Width, Level); ++i) {
            auto j = permute_idx(i);
            if (j > i) {
                std::swap(tile[i], tile[j]);
            }
        }
    }

    constexpr size_type operator()(size_type idx) const noexcept
    {
        return permute_idx(idx);
    }
};

template<unsigned Width>
struct reversing_permutation<Width, 1> {
    using size_type = size_t;
    static const unsigned Level = 1;

    static constexpr size_type permute_idx(size_type idx)
    {
        return idx;
    }

    /// Operate inplace on a single tile
    template<typename T>
    void operator()(T* __restrict tile) const noexcept
    {
        // Do Nothing!
    }

    constexpr size_type operator()(size_type idx) const noexcept
    {
        return permute_idx(idx);
    }
};

template<unsigned Width>
struct reversing_permutation<Width, 0> {
    using size_type = size_t;
    static const unsigned Level = 0;

    static constexpr size_type permute_idx(size_type idx)
    {
        return idx;
    }

    /// Operate inplace on a single tile
    template<typename T>
    void operator()(T* __restrict tile) const noexcept
    {
    }

    constexpr size_type operator()(size_type idx) const noexcept
    {
        return permute_idx(idx);
    }
};



} // namespace dtl
} // namespace alg

#endif//LIBALGEBRA_LIBALGEBRA_DETAIL_REVERSING_PERMUTATION_H_
