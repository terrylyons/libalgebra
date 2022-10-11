/* *************************************************************

Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurkó and Arend Janssen.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */

//  tensor.h

// Include once wrapper
#ifndef DJC_COROPA_LIBALGEBRA_TENSORH_SEEN
#define DJC_COROPA_LIBALGEBRA_TENSORH_SEEN

#ifndef LIBALGEBRA_L1_CACHE_SIZE
#define LIBALGEBRA_L1_CACHE_SIZE 32768// 32Kb should be fairly standard
#endif

//#include <omp.h>

#include <algorithm>
#include <unordered_set>

#include <boost/call_traits.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/functional/hash.hpp>

#include "algebra.h"
#include "area_tensor_basis.h"
#include "base_vector.h"
#include "dense_storage.h"
#include "dense_vector.h"
#include "detail/integer_maths.h"
#include "half_shuffle_tensor_basis.h"
#include "tensor_basis.h"

#define LA_RESTRICT __restrict
#define LA_INLINE_ALWAYS __attribute__((always_inline))

namespace alg {

namespace dtl {

template<unsigned Width, unsigned Level>
struct reversing_permutation {
    using size_type = size_t;
    using cycle_type = std::pair<size_type, size_type>;

    //static const std::vector<cycle_type> cycles;
    static constexpr size_type factor = power(Width, Level - 1);

    static std::vector<cycle_type> make_permutation()
    {
        constexpr size_type tile_size = power(Width, Level);
        std::vector<cycle_type> result;

        std::array<size_type, Level> word{};
        std::array<size_type, Level> rword{};

        auto idx = [](std::array<size_type, Level>& w) {
            size_type result = 0;
            for (auto& v : w) {
                result *= Width;
                result += v;
            }
            return result;
        };

        result.reserve(tile_size);// over-sized, but that's fine
        std::unordered_set<size_type> seen;

        for (size_type lhs = 0; lhs < tile_size; ++lhs) {

            size_type rhs = idx(rword);
            if (lhs != rhs && seen.count(lhs) == 0 && seen.count(rhs) == 0) {
                result.emplace_back(lhs, rhs);
            }

            seen.insert(lhs);
            seen.insert(rhs);

            for (int i = 0; i < Level; ++i) {
                if (word[Level - i - 1] < Width - 1) {
                    ++word[Level - i - 1];
                    ++rword[i];
                    break;
                }
                else {
                    word[Level - i - 1] = 0;
                    rword[i] = 0;
                }
            }
        }

        result.shrink_to_fit();
        return result;
    }

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

    static size_type permute_idx(size_type idx)
    {
        static_assert(Level - 2 > 0, "Level must be at least 3 in this specialisation");
        using next = reversing_permutation<Width, Level - 2>;

        constexpr size_type shift = power(Width, Level - 1);
        return last_letter(idx) * shift + next::permute_idx(middle_word(idx)) * Width + first_letter(idx);
    }

    template<typename T>
    void operator()(T* __restrict tile) const noexcept
    {
        T tmp;

        for (size_type i = 0; i < power(Width, Level); ++i) {
            auto j = permute_idx(i);
            if (j > i) {
                tmp = tile[i];
                tile[i] = tile[j];
                tile[j] = tmp;
            }
        }
    }

    /*
            /// Operate inplace on a single tile
            template <typename T>
            void operator()(T* tile) const noexcept
            {
                T tmp;
                for (auto& cycle : cycles) {
                    tmp = tile[cycle.first];
                    tile[cycle.first] = tile[cycle.second];
                    tile[cycle.second] = tmp;
                }
            }

            /// Operate on different input/ouptut
            template <typename T>
            void operator()(const T* src, T* dst) const noexcept
            {
                for (auto& cycle : cycles) {
                    dst[cycle.first] = src[cycle.second];
                    dst[cycle.second] = src[cycle.first];
                }
            }


            size_type operator()(size_type idx) const
            {
                for (auto& cycle : cycles) {
                    if (idx == cycle.first) {
                        return cycle.second;
                    } else if (idx == cycle.second) {
                        return cycle.first;
                    }
                }
                return idx;
            }
        */
};

template<unsigned Width>
struct reversing_permutation<Width, 2> {
    using size_type = size_t;
    using cycle_type = std::pair<size_type, size_type>;
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
        T tmp;

        for (size_type i = 0; i < power(Width, Level); ++i) {
            auto j = permute_idx(i);
            if (j > i) {
                tmp = tile[i];
                tile[i] = tile[j];
                tile[j] = tmp;
            }
            // tmp = tile[i];
            // tile[i] = tile[j];
            // tile[j] = tmp;
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
    using cycle_type = std::pair<size_type, size_type>;
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
    using cycle_type = std::pair<size_type, size_type>;
    static const unsigned Level = 0;

    static constexpr size_type permute_idx(size_type idx)
    {
        return idx;
    }

    /// Operate inplace on a single tile
    template<typename T>
    void operator()(T* __restrict tile) const noexcept
    {
        // Do nothing!
    }

    constexpr size_type operator()(size_type idx) const noexcept
    {
        return permute_idx(idx);
    }
};

struct default_signer {
    default_signer(DEG degree) : sign(degree % 2 == 0 ? 1 : -1)
    {}

    template<typename S>
    S operator()(const S& arg)
    {
        if (sign == 1) {
            return arg;
        }
        else {
            return -arg;
        }
    }

private:
    int sign;
};

struct non_signing_signer {
    template<typename S>
    S operator()(const S& arg)
    {
        return arg;
    }

    non_signing_signer(DEG deg) {}
};

template<DEG Width, DEG MaxDepth, DEG BlockLetters, typename Scalar, typename Signer>
class tiled_inverse_operator
{
    typedef Scalar SCA;

    using BASIS = free_tensor_basis<Width, MaxDepth>;

public:
    static constexpr DEG block_letters = BlockLetters;
    static constexpr size_t block_width = power(Width, BlockLetters);
    static constexpr size_t block_size = power(Width, 2 * BlockLetters);
    static constexpr unsigned max_middle_word_length = MaxDepth - 2 * BlockLetters;
    //        static constexpr size_type middle_word_count = tensor_alg_size(max_middle_word_length);
    static constexpr size_t block_offset = power(Width, BlockLetters);

    template<DEG Level, DEG MaxLevel>
    struct recursive_untiled_compute {
        using permutation_t = reversing_permutation<Width, Level>;
        //            using signer_t = signing_operator<Level % 2>;
        using next_t = recursive_untiled_compute<Level + 1, MaxLevel>;
        static constexpr size_t level_size = power(Width, Level);

        /*
             * Everything here is supposed to fit in cache, so we really don't need to worry
             * about locality etc when manipulating arrays. This should be very fast.
             */
        void operator()(const SCA* __restrict src_ptr, SCA* __restrict dst_ptr, DEG curr_deg) const noexcept
        {
            assert(curr_deg >= Level);
            // Copy from src to test and adjust sign.
            //                signer_t signer;
            Signer signer(Level);
            for (size_t i = 0; i < level_size; ++i) {
                dst_ptr[i] = signer(src_ptr[i]);
            }

            // Operate on the pointer as if it were a tile of size Width^Level
            permutation_t permutation;
            permutation(dst_ptr);

            // Recurse down to the next level.
            if (curr_deg > Level) {
                next_t next;
                next(src_ptr + level_size, dst_ptr + level_size, curr_deg);
            }
        }
    };

    template<DEG MaxLevel>
    struct recursive_untiled_compute<MaxLevel, MaxLevel> {
        // So the code doesn't change too much, set Level = MaxLevel
        static constexpr DEG Level = MaxLevel;
        using permutation_t = reversing_permutation<Width, Level>;
        //            using signer_t = signing_operator<Level % 2>;
        static constexpr size_t level_size = power(Width, Level);

        /*
             * Everything here is supposed to fit in cache, so we really don't need to worry
             * about locality etc when manipulating arrays. This should be very fast.
             */
        void operator()(const SCA* __restrict src_ptr, SCA* __restrict dst_ptr, DEG) const noexcept
        {
            // Copy from src to test and adjust sign.
            //                signer_t signer;
            Signer signer(Level);
            for (size_t i = 0; i < level_size; ++i) {
                dst_ptr[i] = signer(src_ptr[i]);
            }

            // Operate on the pointer as if it were a tile of size Width^Level
            permutation_t permutation;
            permutation(dst_ptr);

            // No more recursing!
        }
    };

    recursive_untiled_compute<0U, 2 * BlockLetters - 1> recurse;

    void operator()(const SCA* src_ptr, SCA* dst_ptr, const DEG curr_degree) const noexcept
    {

        if (src_ptr == nullptr)// if pointer to source is null
        {
            return;
        }
        if (dst_ptr == nullptr) {
            return;
        }

        recurse(src_ptr, dst_ptr, curr_degree);

        for (unsigned int length = 0; length <= max_middle_word_length && length + 2 * block_letters <= curr_degree; ++length) {
            auto istart = BASIS::start_of_degree(length);
            auto iend = BASIS::start_of_degree(length + 1);

            Signer signer(length);
            auto src_dst_offset = BASIS::start_of_degree(length + 2 * BlockLetters);

            auto src_p = src_ptr + src_dst_offset;
            auto dst_p = dst_ptr + src_dst_offset;

            auto key_start = BASIS::index_to_key(istart);
            auto key_end = BASIS::index_to_key(iend);

            auto word_idx = istart;

            //#pragma omp parallel for
            for (auto word = key_start; word != key_end; word = BASIS::nextkey(word), ++word_idx) {
                auto rword_index = BASIS::key_to_index(word.reverse());

                if (length % 2 == 0) {
                    process_tile(src_p, dst_p, word_idx - istart, rword_index - istart, length, signer);
                }
                else {
                    process_tile(src_p, dst_p, word_idx - istart, rword_index - istart, length, signer);
                }
            }
        }
    }

private:
    static void read_tile(const SCA* __restrict data_ptr, SCA* __restrict tile_ptr, DIMN stride)
    {
        for (size_t row = 0; row < block_width; ++row) {
            auto row_offset = row * stride;
            for (size_t col = 0; col < block_width; ++col) {
                *(tile_ptr++) = data_ptr[row_offset + col];
            }
        }
    }

    static void sign_tile(SCA tile[block_size], Signer signer) noexcept
    {
        for (size_t i = 0; i < block_size; ++i) {
            tile[i] = signer(tile[i]);// tile[i] = op(tile[i]);
        }
    }

    static void write_tile(SCA* __restrict tile_ptr, SCA* __restrict data_ptr, DIMN stride)
    {
        for (size_t row = 0; row < block_width; ++row) {
            auto row_offset = row * stride;

            for (size_t col = 0; col < block_width; ++col) {
                data_ptr[row_offset + col] = *(tile_ptr++);
            }
        }
    }

    void process_tile(
            const SCA* input_data,
            SCA* output_data,
            size_t word_index,
            size_t rword_index,
            DEG degree,
            Signer signer) const
    {
        SCA tile[block_size];
        auto stride = power(Width, degree + BlockLetters);
        assert((block_width - 1) * stride + word_index * block_offset + (block_width - 1) < power(Width, degree + 2 * BlockLetters));
        read_tile(input_data + word_index * block_offset, tile, stride);
        reversing_permutation<Width, 2 * BlockLetters> permutation;
        permutation(tile);
        sign_tile(tile, signer);
        write_tile(tile, output_data + rword_index * block_offset, stride);
    }
};

template<DEG Width, DEG Depth, typename Coeffs>
class free_tensor_multiplication_helper
{

protected:
    using basis_type = tensor_basis<Width, Depth>;
    using coefficient_ring = Coeffs;
    using scalar_type = typename Coeffs::S;
    using pointer = scalar_type*;
    using const_pointer = const scalar_type*;
    using reference = scalar_type&;
    using const_reference = const scalar_type&;

    using tsi = tensor_size_info<Width>;

    template<typename B>
    using dense_tensor = vectors::dense_vector<B, coefficient_ring>;

    const_pointer lhs_data;
    const_pointer rhs_data;
    pointer out_data;

    IDEG lhs_deg, rhs_deg, out_deg, old_out_deg;

    template<typename B>
    void resize_result(dense_tensor<B>& result) const
    {
        // Special treatment for old_out_deg == 0 because otherwise an empty
        // vector is never resized
        if (out_deg > old_out_deg || old_out_deg == 0) {
            result.resize_to_degree(out_deg);
        }
    }

public:
    template<typename B1, typename B2, typename B3>
    free_tensor_multiplication_helper(dense_tensor<B1>& out,
                                      const dense_tensor<B2>& left,
                                      const dense_tensor<B3>& right,
                                      DEG max_degree)
        : lhs_deg(left.degree()), rhs_deg(right.degree()), out_deg(IDEG(max_degree)),
          old_out_deg(IDEG(out.degree()))
    {
        static_assert(std::is_base_of<basis_type, B1>::value
                              && std::is_base_of<basis_type, B2>::value
                              && std::is_base_of<basis_type, B3>::value,
                      "bases are not tensor bases");

        out_deg = std::min(out_deg, lhs_deg + rhs_deg);
        resize_result(out);

        lhs_data = left.as_ptr();
        rhs_data = right.as_ptr();
        out_data = out.as_mut_ptr();
    }

    template<typename B1, typename B2>
    free_tensor_multiplication_helper(dense_tensor<B1>& left,
                                      const dense_tensor<B2>& right,
                                      DEG max_degree)
        : lhs_deg(left.degree()), rhs_deg(right.degree()), out_deg(IDEG(max_degree)),
          old_out_deg(IDEG(lhs_deg))
    {
        out_deg = std::min(out_deg, lhs_deg + rhs_deg);
        resize_result(left);

        lhs_data = nullptr;
        rhs_data = right.as_ptr();
        out_data = left.as_mut_ptr();
    }

    IDEG old_out_degree() const noexcept { return old_out_deg; }
    IDEG lhs_degree() const noexcept { return lhs_deg; }
    IDEG rhs_degree() const noexcept { return rhs_deg; }
    IDEG out_degree() const noexcept { return out_deg; }

    reference out_unit() noexcept { return out_data[0]; }
    const_reference left_unit() const noexcept { return (lhs_data == nullptr) ? out_data[0] : lhs_data[0]; }
    const_reference right_unit() const noexcept { return rhs_data[0]; }
    const_pointer left_fwd_read(IDEG d, IDIMN offset = 0) const noexcept
    {
        assert(d >= 0 && d <= lhs_deg);
        assert(offset + basis_type::start_of_degree(d) <= basis_type::start_of_degree(d + 1));
        return ((lhs_data == nullptr) ? out_data : lhs_data) + offset + basis_type::start_of_degree(d);
    }
    const_pointer right_fwd_read(IDEG d, IDIMN offset = 0) const noexcept
    {
        assert(d >= 0 && d <= rhs_deg);
        assert(offset + basis_type::start_of_degree(d) <= basis_type::start_of_degree(d + 1));
        return rhs_data + offset + basis_type::start_of_degree(d);
    }
    pointer fwd_write(IDEG d) const noexcept
    {
        assert(d >= 0 && d <= out_deg);
        return out_data + basis_type::start_of_degree(d);
    }

    std::pair<DIMN, DIMN> range_size(IDEG lhs, IDEG rhs) const noexcept
    {
        return std::pair<DIMN, DIMN>(tsi::powers[lhs], tsi::powers[rhs]);
    }
};

template<DEG Width, DEG Depth, IDEG TileLetters = 0>
struct tile_details {
    static constexpr IDEG tile_letters = TileLetters;
    static constexpr DIMN tile_width = integer_maths::power(Width, tile_letters);
    static constexpr DIMN tile_size = tile_width * tile_width;
    static constexpr DIMN tile_shift = integer_maths::power(Width, tile_letters - 1);
};

template<DEG Width, DEG Depth, typename Coeffs, IDEG TileLetters = 1>
class tiled_free_tensor_multiplication_helper
    : public free_tensor_multiplication_helper<Width, Depth, Coeffs>,
      public tile_details<Width, Depth, TileLetters>
{
    using base = free_tensor_multiplication_helper<Width, Depth, Coeffs>;
    using tile_info = tile_details<Width, Depth, TileLetters>;
    using tsi = tensor_size_info<Width>;

public:
    using scalar_type = typename Coeffs::SCA;
    using pointer = scalar_type*;
    using const_pointer = const scalar_type*;
    using reference = scalar_type&;
    using const_reference = const scalar_type&;

private:
    template<typename B>
    using dense_tensor = vectors::dense_vector<B, Coeffs>;

    using basis_type = tensor_basis<Width, Depth>;

    std::vector<scalar_type> left_read_tile;
    std::vector<scalar_type> right_read_tile;
    std::vector<scalar_type> output_tile;
    std::vector<scalar_type> reverse_data;
    const_pointer left_reverse_read_ptr = nullptr;
    pointer reverse_write_ptr = nullptr;

public:
    using tile_info::tile_letters;
    using tile_info::tile_shift;
    using tile_info::tile_size;
    using tile_info::tile_width;

    template<typename B1, typename B2, typename B3>
    tiled_free_tensor_multiplication_helper(dense_tensor<B1>& out,
                                            const dense_tensor<B2>& lhs,
                                            const dense_tensor<B3>& rhs,
                                            DEG max_degree)
        : base(out, lhs, rhs, max_degree),
          left_read_tile(tile_width),
          right_read_tile(tile_width),
          output_tile(tile_size)
    {
        if (base::lhs_deg > 0) {
            reverse_data.resize(basis_type::start_of_degree(base::lhs_deg + 1));
            dtl::tiled_inverse_operator<Width, (Depth > 0) ? Depth - 1 : 0, TileLetters, scalar_type, dtl::non_signing_signer> reverser;
            reverser(base::lhs_data, reverse_data.data(), base::lhs_deg);
            left_reverse_read_ptr = reverse_data.data();
        }
    }

    template<typename B1, typename B2>
    tiled_free_tensor_multiplication_helper(dense_tensor<B1>& lhs, const dense_tensor<B2>& rhs, DEG max_degree)
        : base(lhs, rhs, max_degree),
          left_read_tile(tile_width),
          right_read_tile(tile_width),
          output_tile(tile_size)
    {
        if (base::lhs_deg > 0) {
            reverse_data.resize(basis_type::start_of_degree(base::lhs_deg + 1));
            dtl::tiled_inverse_operator<Width, (Depth > 0) ? Depth - 1 : 0, TileLetters, scalar_type, dtl::non_signing_signer> reverser;
            reverser(base::out_data, reverse_data.data(), base::lhs_deg);
            left_reverse_read_ptr = reverse_data.data();
        }
    }

/*    template<typename B>
    tiled_free_tensor_multiplication_helper(
            dense_tensor<free_tensor_basis<Width, Depth>>& out,
            const dense_tensor<free_tensor_basis<Width, Depth>>& lhs,
            const dense_tensor<B>& rhs,
            DEG max_degree)
        : base(out, lhs, rhs, max_degree),
          left_read_tile(tile_width),
          right_read_tile(tile_width),
          output_tile(tile_size)
    {
//        assert(base::out_deg == 0 || out.reverse_dimension() == tsi::degree_sizes[base::out_deg-1]);
//        reverse_write_ptr = out.as_mut_rptr();
//        left_reverse_read_ptr = lhs.as_rptr();


        if (base::lhs_deg > 0) {
            reverse_data.resize(basis_type::start_of_degree(base::lhs_deg + 1));
            dtl::tiled_inverse_operator<Width, (Depth > 0) ? Depth - 1 : 0, TileLetters, scalar_type, dtl::non_signing_signer> reverser;
            reverser(base::out_data, reverse_data.data(), base::lhs_deg);
            left_reverse_read_ptr = reverse_data.data();
        }
    }

    template<typename B>
    tiled_free_tensor_multiplication_helper(
            dense_tensor<free_tensor_basis<Width, Depth>>& lhs,
            const dense_tensor<B>& rhs,
            DEG max_degree)
        : base(lhs, rhs, max_degree),
          left_read_tile(tile_width),
          right_read_tile(tile_width),
          output_tile(tile_size)
    {
//        assert(base::out_deg == 0 || lhs.reverse_dimension() >= tsi::degree_sizes[base::out_deg-1]);
//        reverse_write_ptr = lhs.as_mut_rptr();
//        left_reverse_read_ptr = reverse_write_ptr;
        if (base::lhs_deg > 0) {
            reverse_data.resize(basis_type::start_of_degree(base::lhs_deg + 1));
            dtl::tiled_inverse_operator<Width, (Depth > 0) ? Depth - 1 : 0, TileLetters, scalar_type, dtl::non_signing_signer> reverser;
            reverser(base::out_data, reverse_data.data(), base::lhs_deg);
            left_reverse_read_ptr = reverse_data.data();
        }
    }*/

    pointer out_tile_ptr() noexcept
    {
        return output_tile.data();
    }
    const_pointer left_read_tile_ptr() const noexcept
    {
        return left_read_tile.data();
    }
    const_pointer right_read_tile_ptr() const noexcept
    {
        return right_read_tile.data();
    }

    void read_left_tile(IDEG degree, IDIMN index) noexcept
    {
        const auto start_of_degree = basis_type::start_of_degree(degree);
        const auto* ptr_begin = left_reverse_read_ptr + index * tile_width + start_of_degree;
        std::copy(ptr_begin, ptr_begin + tile_width, left_read_tile.data());
    }
    void read_right_tile(IDEG degree, IDIMN index) noexcept
    {
        const auto start_of_degree = basis_type::start_of_degree(degree);
        const auto* ptr_begin = base::rhs_data + index * tile_width + start_of_degree;
        std::copy(ptr_begin, ptr_begin + tile_width, right_read_tile.data());
    }

    void reset_tile(IDEG degree, IDIMN index, IDIMN /*reverse_index*/) noexcept
    {
        assert(0 <= degree && degree <= static_cast<IDEG>(Depth));
        assert(index < static_cast<IDEG>(tsi::powers[degree - 2 * TileLetters]));
        const auto start_of_degree = basis_type::start_of_degree(degree);
        const auto stride = tsi::powers[degree - tile_letters];

        const_pointer optr = base::out_data + index * tile_width + start_of_degree;
        pointer tptr = output_tile.data();

        for (DIMN i = 0; i < tile_width; ++i) {
            for (DIMN j = 0; j < tile_width; ++j) {
                tptr[i * tile_width + j] = optr[i * stride + j];
            }
        }
    }

    void reset_tile_to_zero() noexcept
    {
        std::fill(output_tile.begin(), output_tile.end(), Coeffs::zero);
    }

    void write_tile(IDEG degree, IDIMN index, IDIMN reverse_index) noexcept
    {
        assert(0 <= degree && degree <= static_cast<IDEG>(Depth));
        assert(index <= static_cast<IDIMN>(tsi::powers[degree - 2 * TileLetters]));
        assert(reverse_index <= static_cast<IDIMN>(tsi::powers[degree - 2 * TileLetters]));
        const auto start_of_degree = basis_type::start_of_degree(degree);
        pointer optr = base::out_data + index * tile_width + start_of_degree;
        const_pointer tptr = output_tile.data();
        auto stride = tsi::powers[degree - tile_letters];

        for (DIMN i = 0; i < tile_width; ++i) {
            for (DIMN j = 0; j < tile_width; ++j) {
                optr[i * stride + j] = tptr[i * tile_width + j];
            }
        }

        if (reverse_write_ptr != nullptr && degree < base::out_deg) {
            // Write out reverse data
//            using perm = reversing_permutation<Width, tile_info::tile_letters>;
//
//            assert(((tile_width-1)*stride + (tile_width-1) + reverse_index*tile_width+start_of_degree) < tsi::degree_sizes[degree]);
//            optr = reverse_write_ptr + reverse_index * tile_width + start_of_degree;
//            for (DIMN i = 0; i < tile_width; ++i) {
//                for (DIMN j = 0; j < tile_width; ++j) {
//                    optr[i * stride + j] = tptr[perm::permute_idx(i) * tile_width + j];
//                }
//            }
        }
    }

    static std::pair<IDIMN, IDIMN> split_key(IDEG split_deg, IDIMN key) noexcept
    {
        const auto splitter = static_cast<IDIMN>(tsi::powers[split_deg]);
        return {key / splitter, key % splitter};
    }

    static IDIMN combine_keys(IDEG right_deg, IDIMN left, IDIMN right) noexcept
    {
        const auto shift = static_cast<IDIMN>(tsi::powers[right_deg]);
        return left * shift + right;
    }

    static IDIMN reverse_key(IDEG degree, IDIMN index) noexcept
    {
        IDIMN result = 0;
        for (auto i = 0; i < degree; ++i) {
            result *= Width;
            result += (index % Width);
            index /= Width;
        }
        return result;
        //        if (degree < 2) {
        //            return index;
        //        }
        //        if (degree == 2) {
        //            return (index % Width) * Width + (index / Width);
        //        }
        //
        //        assert(degree < static_cast<IDEG>(tsi::powers.size()));
        //        const auto high = static_cast<IDIMN>(tsi::powers[degree-1]);
        //        auto left_letter = index / high;
        //        auto right_letter = index % Width;
        //        auto middle_letters = (index % high) / Width;
        //        return combine_keys(1,
        //                            combine_keys(degree-1,
        //                                         right_letter,
        //                                         reverse_key(degree - 2, middle_letters)),
        //                            left_letter);
    }

    const_pointer left_fwd_read_ptr(IDEG degree, IDIMN index) const noexcept
    {
        return base::left_fwd_read(degree, index * tile_width);
    }
    const_pointer right_fwd_read_ptr(IDEG degree, IDIMN index) const noexcept
    {
        return base::right_fwd_read(degree, index * tile_width);
    }
};

}// namespace dtl

template<DEG Width, DEG Depth>
class free_tensor_multiplier
{
public:
    using basis_type = tensor_basis<Width, Depth>;
    using key_type = typename basis_type::KEY;
    using pair_type = std::pair<key_type, int>;
    using result_type = boost::container::small_vector<pair_type, 1>;
    using argument_type = key_type;

    result_type operator()(argument_type lhs, argument_type rhs) const
    {
        assert(lhs.valid() && rhs.valid());
        return {{lhs * rhs, 1}};
    }
};

template<DEG Width, DEG Depth>
class traditional_free_tensor_multiplication
    : public base_multiplication<free_tensor_multiplier<Width, Depth>,
                                 traditional_free_tensor_multiplication<Width, Depth>>,
      public dtl::hybrid_vector_mixin<traditional_free_tensor_multiplication<Width, Depth>>
{
    using base = base_multiplication<free_tensor_multiplier<Width, Depth>,
                                     traditional_free_tensor_multiplication>;
    using mixin = dtl::hybrid_vector_mixin<traditional_free_tensor_multiplication>;

    using tsi = dtl::tensor_size_info<Width>;

protected:
    using base::m_multiplier;

public:
    using base::fma;
    using base::multiply_inplace;
    using mixin::fma;
    using mixin::multiply_inplace;

    using basis_type = tensor_basis<Width, Depth>;

protected:
    template<typename B, typename Coeffs>
    static void update_reverse_data(vectors::dense_vector<B, Coeffs>& out, DEG max_degree)
    {}

    template<typename Coeffs>
    static void update_reverse_data(vectors::dense_vector<free_tensor_basis<Width, Depth>, Coeffs>& out, DEG max_degree)
    {
        if (max_degree > 0) {
            out.construct_reverse_data(max_degree-1);
        }
    }

private:
    template<typename Coeffs, typename Fn>
    LA_INLINE_ALWAYS static void
    update_inplace_lhs(typename Coeffs::S* LA_RESTRICT lhs_ptr,
                       const typename Coeffs::S& rhs_val,
                       IDIMN lhs_count,
                       Fn op) noexcept
    {
        for (IDIMN i = 0; i < lhs_count; ++i) {
            lhs_ptr[i] = op(Coeffs::template mul(lhs_ptr[i], rhs_val));
        }
    }

    template<typename Coeffs, typename Fn>
    LA_INLINE_ALWAYS static void update_rhs_max(typename Coeffs::S* LA_RESTRICT lhs_ptr,
                                                const typename Coeffs::S& lhs_unit,
                                                const typename Coeffs::S* LA_RESTRICT rhs_ptr,
                                                IDIMN rhs_count,
                                                Fn op) noexcept
    {
        for (IDIMN i = 0; i < rhs_count; ++i) {
            Coeffs::template add_inplace(lhs_ptr[i], op(Coeffs::template mul(lhs_unit, rhs_ptr[i])));
        }
    }

    template<typename Coeffs, typename Fn>
    LA_INLINE_ALWAYS static void assign_rhs_max(typename Coeffs::S* LA_RESTRICT lhs_ptr,
                                                const typename Coeffs::S& lhs_unit,
                                                const typename Coeffs::S* LA_RESTRICT rhs_ptr,
                                                IDIMN rhs_count,
                                                Fn op) noexcept
    {
        for (IDIMN i = 0; i < rhs_count; ++i) {
            lhs_ptr[i] = op(Coeffs::template mul(lhs_unit, rhs_ptr[i]));
        }
    }

    template<typename Coeffs, typename Fn>
    LA_INLINE_ALWAYS static void update_accumulate(typename Coeffs::S* LA_RESTRICT optr,
                                                   const typename Coeffs::S* LA_RESTRICT lhs_ptr,
                                                   const typename Coeffs::S* LA_RESTRICT rhs_ptr,
                                                   IDIMN lhs_count,
                                                   IDIMN rhs_count,
                                                   Fn op) noexcept
    {
        for (IDIMN i = 0; i < lhs_count; ++i) {
            for (IDIMN j = 0; j < rhs_count; ++j) {
                Coeffs::template add_inplace(*(optr++), op(Coeffs::template mul(lhs_ptr[i], rhs_ptr[j])));
            }
        }
    }

    template<typename Coeffs, typename Fn>
    LA_INLINE_ALWAYS static void
    update_assign(typename Coeffs::S* LA_RESTRICT optr,
                  const typename Coeffs::S* LA_RESTRICT lhs_ptr,
                  const typename Coeffs::S* LA_RESTRICT rhs_ptr,
                  IDIMN lhs_count,
                  IDIMN rhs_count,
                  Fn op) noexcept
    {
        for (IDIMN i = 0; i < lhs_count; ++i) {
            for (IDIMN j = 0; j < rhs_count; ++j) {
                *(optr++) = op(Coeffs::template mul(lhs_ptr[i], rhs_ptr[j]));
            }
        }
    }

protected:
    template<typename Coeffs>
    using helper = dtl::free_tensor_multiplication_helper<Width, Depth, Coeffs>;

    template<typename Coeffs, typename Fn>
    void fma_impl(helper<Coeffs>& helper, Fn op, IDEG max_degree) const noexcept
    {
        for (IDEG out_deg = max_degree; out_deg > 0; --out_deg) {
            auto lhs_deg_min = std::max(0, out_deg - helper.rhs_degree());
            auto lhs_deg_max = std::min(out_deg, helper.lhs_degree());

            auto* out_ptr = helper.fwd_write(out_deg);

            for (auto lh_deg = lhs_deg_max; lh_deg >= lhs_deg_min; --lh_deg) {
                auto rh_deg = out_deg - lh_deg;
                auto lhs_ptr = helper.left_fwd_read(lh_deg);
                auto rhs_ptr = helper.right_fwd_read(rh_deg);

                auto* p = out_ptr;
                auto sizes = helper.range_size(lh_deg, rh_deg);

                for (IDIMN i = 0; i < IDIMN(sizes.first); ++i) {
                    for (IDIMN j = 0; j < IDIMN(sizes.second); ++j) {
                        Coeffs::template add_inplace(*(p++),
                                                     op(Coeffs::template mul(lhs_ptr[i], rhs_ptr[j])));
                    }
                }
            }
        }

        auto& out_unit = helper.out_unit();
        out_unit = op(Coeffs::template mul(helper.left_unit(), helper.right_unit()));
    }

    template<typename Coeffs, typename Fn>
    void multiply_inplace_impl(helper<Coeffs>& helper, Fn op, IDEG max_degree) const noexcept
    {

        const auto old_out_deg = helper.lhs_degree();
        const auto rhs_max_deg = helper.rhs_degree();

        const auto& lhs_unit = helper.left_unit();
        const auto& rhs_unit = helper.right_unit();
        bool default_assign = rhs_unit == Coeffs::zero;

        for (IDEG out_deg = max_degree; out_deg > 0; --out_deg) {
            bool assign = out_deg > old_out_deg || default_assign;
            auto lhs_deg_min = std::max(IDEG(1), out_deg - rhs_max_deg);
            auto lhs_deg_max = std::min(out_deg - 1, old_out_deg);

            if (!assign) {
                update_inplace_lhs<Coeffs>(helper.fwd_write(out_deg),
                                           rhs_unit,
                                           static_cast<IDIMN>(tsi::powers[out_deg]),
                                           op);
            }
            else if (default_assign && out_deg > rhs_max_deg && lhs_deg_max < lhs_deg_min) {
                update_inplace_lhs<Coeffs>(helper.fwd_write(out_deg), rhs_unit,
                                           static_cast<IDIMN>(tsi::powers[out_deg]), op);
            }

            for (IDEG lhs_deg = lhs_deg_max; lhs_deg >= lhs_deg_min; --lhs_deg) {
                auto rhs_deg = out_deg - lhs_deg;

                const auto lhs_count = static_cast<IDIMN>(tsi::powers[lhs_deg]);
                const auto rhs_count = static_cast<IDIMN>(tsi::powers[rhs_deg]);

                if (assign) {
                    update_assign<Coeffs>(helper.fwd_write(out_deg),
                                          helper.left_fwd_read(lhs_deg),
                                          helper.right_fwd_read(rhs_deg),
                                          lhs_count,
                                          rhs_count,
                                          op);
                    assign = false;
                }
                else {
                    update_accumulate<Coeffs>(helper.fwd_write(out_deg),
                                              helper.left_fwd_read(lhs_deg),
                                              helper.right_fwd_read(rhs_deg),
                                              lhs_count,
                                              rhs_count,
                                              op);
                }
            }

            if (out_deg <= rhs_max_deg) {
                if (assign) {
                    assign_rhs_max<Coeffs>(helper.fwd_write(out_deg),
                                           lhs_unit,
                                           helper.right_fwd_read(out_deg),
                                           static_cast<IDIMN>(tsi::powers[out_deg]),
                                           op);
                }
                else {
                    update_rhs_max<Coeffs>(helper.fwd_write(out_deg),
                                           lhs_unit,
                                           helper.right_fwd_read(out_deg),
                                           static_cast<IDIMN>(tsi::powers[out_deg]),
                                           op);
                }
            }
        }

        auto& out_unit = helper.out_unit();
        out_unit = op(Coeffs::template mul(lhs_unit, rhs_unit));
    }

public:
    template<typename Basis, typename Coeffs, typename Fn, typename OriginalVectors>
    std::enable_if_t<std::is_base_of<basis_type, Basis>::value>
    fma(vectors::dense_vector<Basis, Coeffs>& out,
        const vectors::dense_vector<Basis, Coeffs>& lhs,
        const vectors::dense_vector<Basis, Coeffs>& rhs,
        Fn op,
        DEG max_degree,
        OriginalVectors&) const
    {
        if (!lhs.empty() && !rhs.empty()) {
            helper<Coeffs> help(out, lhs, rhs, max_degree);
            fma_impl(help, op, help.out_degree());
//            update_reverse_data(out, help.out_degree());
        }
    }

    template<typename Basis, typename Coeffs, typename Fn, typename OriginalVectors>
    std::enable_if_t<std::is_base_of<basis_type, Basis>::value>
    multiply_inplace(vectors::dense_vector<Basis, Coeffs>& lhs,
                     const vectors::dense_vector<Basis, Coeffs>& rhs,
                     Fn op,
                     DEG max_degree, OriginalVectors&) const
    {
        if (!rhs.empty()) {
            helper<Coeffs> help(lhs, rhs, max_degree);
            multiply_inplace_impl(help, op, help.out_degree());
//            update_reverse_data(lhs, help.out_degree());
        }
        else {
            lhs.clear();
        }
    }
};

template<DEG Width, DEG Depth, DEG TileLetters = 1>
class tiled_free_tensor_multiplication
    : public traditional_free_tensor_multiplication<Width, Depth>
{
    using base = traditional_free_tensor_multiplication<Width, Depth>;

    template<typename C>
    using helper_type = dtl::tiled_free_tensor_multiplication_helper<
            Width, Depth, C, TileLetters>;

    using tile_info = dtl::tile_details<Width, Depth, TileLetters>;
    using tsi = dtl::tensor_size_info<Width>;

    template<typename S, typename Fn>
    LA_INLINE_ALWAYS static void impl_0bd(S* LA_RESTRICT tile,
                                          const S& lhs_unit,
                                          const S* LA_RESTRICT rhs_ptr,
                                          IDIMN stride,
                                          Fn op) noexcept
    {
        constexpr auto tile_width = static_cast<IDIMN>(tile_info::tile_width);
        for (IDIMN i = 0; i < tile_width; ++i) {
            for (IDIMN j = 0; j < tile_width; ++j) {
                tile[i * tile_width + j] += op(lhs_unit * rhs_ptr[i * stride + j]);
            }
        }
    }

    template<typename S, typename Fn>
    LA_INLINE_ALWAYS static void impl_db0(S* LA_RESTRICT tile,
                                          const S* LA_RESTRICT lhs_ptr,
                                          const S& rhs_unit,
                                          IDIMN stride,
                                          Fn op) noexcept
    {
        constexpr auto tile_width = static_cast<IDIMN>(tile_info::tile_width);
        for (IDIMN i = 0; i < tile_width; ++i) {
            for (IDIMN j = 0; j < tile_width; ++j) {
                tile[i * tile_width + j] += op(lhs_ptr[i * stride + j] * rhs_unit);
            }
        }
    }

    template<typename S, typename Fn>
    LA_INLINE_ALWAYS static void impl_mid(S* LA_RESTRICT tile,
                                          const S* LA_RESTRICT lhs_tile,
                                          const S* LA_RESTRICT rhs_tile,
                                          Fn op) noexcept
    {
        using perm = dtl::reversing_permutation<Width, tile_info::tile_letters>;
        constexpr auto tile_width = static_cast<IDIMN>(tile_info::tile_width);
        for (IDIMN i = 0; i < tile_width; ++i) {
            for (IDIMN j = 0; j < tile_width; ++j) {
                tile[perm::permute_idx(i) * tile_width + j] += op(lhs_tile[i] * rhs_tile[j]);
            }
        }
    }

    template<typename S, typename Fn>
    LA_INLINE_ALWAYS static void impl_lb1(S* LA_RESTRICT tile,
                                          const S& lhs_val,
                                          const S* LA_RESTRICT rhs_tile,
                                          Fn op,
                                          IDIMN i) noexcept
    {
        constexpr auto tile_width = static_cast<IDIMN>(tile_info::tile_width);
        for (IDIMN j = 0; j < tile_width; ++j) {
            tile[i * tile_width + j] += op(lhs_val * rhs_tile[j]);
        }
    }

    template<typename S, typename Fn>
    LA_INLINE_ALWAYS static void impl_1br(S* LA_RESTRICT tile,
                                          const S* LA_RESTRICT lhs_tile,
                                          const S& rhs_val,
                                          Fn op,
                                          IDIMN j) noexcept
    {
        using perm = dtl::reversing_permutation<Width, tile_info::tile_letters>;
        constexpr auto tile_width = static_cast<IDIMN>(tile_info::tile_width);
        for (IDIMN i = 0; i < tile_width; ++i) {
            tile[perm::permute_idx(i) * tile_width + j] += op(lhs_tile[i] * rhs_val);
        }
    }

protected:
    template<typename Coeffs, typename Fn>
    LA_INLINE_ALWAYS void impl_common(helper_type<Coeffs>& helper,
                                      Fn op,
                                      IDEG out_deg,
                                      IDIMN k,
                                      IDIMN k_reverse,
                                      IDEG lhs_deg_min,
                                      IDEG lhs_deg_max) const
    {
        constexpr auto tile_letters = static_cast<IDEG>(tile_info::tile_letters);
        constexpr auto tile_width = static_cast<IDIMN>(tile_info::tile_width);

        auto* LA_RESTRICT tile = helper.out_tile_ptr();
        const auto* LA_RESTRICT left_rtile = helper.left_read_tile_ptr();
        const auto* LA_RESTRICT right_rtile = helper.right_read_tile_ptr();

        //        const auto& lhs_unit = helper.left_unit();
        //        const auto& rhs_unit = helper.right_unit();

        const auto mid_deg = out_deg - 2 * tile_letters;

        for (IDEG lhs_deg = lhs_deg_min; lhs_deg < std::min(tile_letters, lhs_deg_max); ++lhs_deg) {
            const auto rhs_deg = out_deg - lhs_deg;

            assert(1 <= lhs_deg && lhs_deg <= tile_letters);
            for (IDIMN i = 0; i < tile_width; ++i) {
                const auto split = helper.split_key(lhs_deg, i);
                const auto& left_val = *helper.left_fwd_read(lhs_deg, split.first);
                helper.read_right_tile(rhs_deg, helper.combine_keys(mid_deg, split.second, k));
                impl_lb1(tile, left_val, right_rtile, op, i);
            }
        }

        for (IDEG lhs_deg = std::max(tile_letters, lhs_deg_min);
             lhs_deg <= std::min(out_deg - tile_letters, lhs_deg_max); ++lhs_deg) {
            const auto rhs_deg = out_deg - lhs_deg;
            assert(tile_letters <= lhs_deg && lhs_deg <= out_deg - tile_letters);
            assert(tile_letters <= rhs_deg && rhs_deg <= out_deg - tile_letters);
            const auto lhs_split = lhs_deg - tile_letters;
            const auto rhs_split = rhs_deg - tile_letters;
            assert(lhs_split + rhs_split == mid_deg);

            const auto split = helper.split_key(rhs_split, k);
            helper.read_left_tile(lhs_deg, helper.reverse_key(lhs_split, split.first));
            helper.read_right_tile(rhs_deg, split.second);
            impl_mid(tile, left_rtile, right_rtile, op);
        }

        for (IDEG lhs_deg = std::max(out_deg - tile_letters + 1, lhs_deg_min);
             lhs_deg <= lhs_deg_max; ++lhs_deg) {
            const auto rhs_deg = out_deg - lhs_deg;
            assert(mid_deg < lhs_deg && lhs_deg < out_deg);
            assert(1 <= rhs_deg && rhs_deg < tile_letters);
            const auto split_left_letters = tile_letters - rhs_deg;
            assert(0 < split_left_letters && split_left_letters < tile_letters);
            assert(lhs_deg == mid_deg + split_left_letters + tile_letters);

            for (IDIMN j = 0; j < tile_width; ++j) {
                const auto split = helper.split_key(rhs_deg, j);
                const auto& right_val = *helper.right_fwd_read(rhs_deg, split.second);
                helper.read_left_tile(lhs_deg, helper.combine_keys(split_left_letters, helper.reverse_key(split_left_letters, split.first), k_reverse));
                impl_1br(tile, left_rtile, right_val, op, j);
            }
        }

        helper.write_tile(out_deg, k, k_reverse);
    }

    template<typename Coeffs, typename Fn>
    void fma_impl(helper_type<Coeffs>& helper, Fn op, DEG max_degree) const
    {
        constexpr auto tile_letters = static_cast<IDEG>(tile_info::tile_letters);

        auto* LA_RESTRICT tile = helper.out_tile_ptr();

        for (IDEG out_deg = static_cast<IDEG>(max_degree); out_deg > 2 * tile_letters; --out_deg) {
            const auto mid_deg = out_deg - 2 * tile_letters;
            const auto stride = static_cast<IDIMN>(tsi::powers[out_deg - tile_letters]);

            auto lhs_deg_min = std::max(IDEG(1), out_deg - helper.rhs_degree());
            auto lhs_deg_max = std::min(out_deg - 1, helper.lhs_degree());

            for (IDIMN k = 0; k < static_cast<IDIMN>(tsi::powers[mid_deg]); ++k) {
                auto k_reverse = helper.reverse_key(mid_deg, k);

                helper.reset_tile(out_deg, k, k_reverse);

                const auto& lhs_unit = helper.left_unit();
                if (helper.rhs_degree() >= out_deg) {
                    impl_0bd(tile, lhs_unit, helper.right_fwd_read_ptr(out_deg, k), stride, op);
                }

                const auto& rhs_unit = helper.right_unit();
                if (helper.lhs_degree() >= out_deg) {
                    impl_db0(tile, helper.left_fwd_read_ptr(out_deg, k), rhs_unit, stride, op);
                }

                impl_common(helper, op, out_deg, k, k_reverse, lhs_deg_min, lhs_deg_max);
            }
        }

        base::fma_impl(helper, op, std::min(max_degree, DEG(2 * tile_info::tile_letters)));
    }

    template<typename Coeffs, typename Fn>
    void multiply_inplace_impl(helper_type<Coeffs>& helper, Fn op, DEG max_degree) const
    {
        constexpr auto tile_letters = static_cast<IDEG>(tile_info::tile_letters);

        auto* LA_RESTRICT tile = helper.out_tile_ptr();

        const auto& lhs_unit = helper.left_unit();
        const auto& rhs_unit = helper.right_unit();

        const auto old_lhs_deg = helper.lhs_degree();
        const auto rhs_max_deg = helper.rhs_degree();

        for (IDEG out_deg = static_cast<IDEG>(max_degree); out_deg > 2 * tile_letters; --out_deg) {
            const auto mid_deg = out_deg - 2 * tile_letters;
            const auto stride = static_cast<IDIMN>(tsi::powers[out_deg - tile_letters]);

            auto lhs_deg_min = std::max(IDEG(1), out_deg - rhs_max_deg);
            auto lhs_deg_max = std::min(out_deg - 1, old_lhs_deg);

            for (IDIMN k = 0; k < static_cast<IDIMN>(tsi::powers[mid_deg]); ++k) {
                auto k_reverse = helper.reverse_key(mid_deg, k);

                //                if (out_deg > old_lhs_deg) {
                //                    helper.reset_tile(out_deg, k, k_reverse);
                //                }
                //                else {
                helper.reset_tile_to_zero();
                //                }
                if (out_deg <= old_lhs_deg && rhs_unit != Coeffs::zero) {
                    impl_db0(tile, helper.left_fwd_read_ptr(out_deg, k), rhs_unit, stride, op);
                }

                if (out_deg <= rhs_max_deg && lhs_unit != Coeffs::zero) {
                    impl_0bd(tile, lhs_unit, helper.right_fwd_read_ptr(out_deg, k), stride, op);
                }
                impl_common(helper, op, out_deg, k, k_reverse, lhs_deg_min, lhs_deg_max);
            }
        }

        base::multiply_inplace_impl(helper, op, std::min(IDEG(max_degree), 2 * tile_info::tile_letters));
    }

public:
    using base::fma;
    using base::multiply_inplace;

#pragma clang diagnostic push
#pragma ide diagnostic ignored "HidingNonVirtualFunction"
    template<typename B1, typename Coeffs, typename Fn, typename OriginalVectors>
    void fma(vectors::dense_vector<B1, Coeffs>& out,
             const vectors::dense_vector<B1, Coeffs>& lhs,
             const vectors::dense_vector<B1, Coeffs>& rhs,
             Fn op,
             DEG max_degree,
             OriginalVectors& orig) const
    {
        if (max_degree <= 2 * tile_info::tile_letters) {
            base::fma(out, lhs, rhs, op, max_degree, orig);
            return;
        }

        if (!lhs.empty() && !rhs.empty()) {
            helper_type<Coeffs> helper(out, lhs, rhs, max_degree);
            fma_impl(helper, op, helper.out_degree());
            if (helper.out_degree() > 0) {
                auto update_deg = std::min(helper.out_degree()-1, 2*tile_info::tile_letters);
//                base::update_reverse_data(out, update_deg);
            }
        }
    }

    template<typename Basis, typename Coeffs, typename Fn, typename OriginalVectors>
    void multiply_inplace(vectors::dense_vector<Basis, Coeffs>& lhs,
                          const vectors::dense_vector<Basis, Coeffs>& rhs,
                          Fn op, DEG max_degree, OriginalVectors& orig) const
    {
        //        std::cout << "BEFORE " << lhs << '\n';
        if (max_degree <= 2 * tile_info::tile_letters) {
            base::multiply_inplace(lhs, rhs, op, max_degree, orig);
            return;
        }

        if (!rhs.empty()) {
            helper_type<Coeffs> helper(lhs, rhs, max_degree);
            multiply_inplace_impl(helper, op, helper.out_degree());
            if (helper.out_degree() > 0) {
                auto update_deg = std::min(helper.out_degree()-1, 2*tile_info::tile_letters);
//                base::update_reverse_data(lhs, update_deg);
            }
        }
        else {
            lhs.clear();
        }
    }
#pragma clang diagnostic pop
};

//template <DEG Width, DEG Depth>
//using free_tensor_multiplication = traditional_free_tensor_multiplication<Width, Depth>;

template<DEG Width, DEG Depth>
using free_tensor_multiplication = tiled_free_tensor_multiplication<Width, Depth, 1>;

template<DEG Width, DEG Depth>
class left_half_shuffle_multiplier
    : public multiplier_base<left_half_shuffle_multiplier<Width, Depth>,
                             tensor_basis<Width, Depth>>
{
    using base = multiplier_base<left_half_shuffle_multiplier, tensor_basis<Width, Depth>>;
    friend base;

    free_tensor_multiplier<Width, Depth> ftmul;

public:
    using basis_type = tensor_basis<Width, Depth>;
    using key_type = typename basis_type::KEY;
    using typename base::argument_type;
    using typename base::inner_result_type;
    using typename base::pair_type;
    using typename base::result_type;

protected:
    inner_result_type shuffle(argument_type lhs, argument_type rhs) const
    {
        return base::add(operator()(lhs, rhs), operator()(rhs, lhs));
    }

private:
    inner_result_type prod_impl(argument_type lhs, argument_type rhs) const
    {
        assert(lhs.valid() && rhs.valid());
        if (rhs.size() == 0) {
            return {{key_type(lhs), 1}};
        }

        result_type recursed = shuffle(lhs.rparent(), rhs);
        const auto k1l(lhs.lparent());

        inner_result_type result;
        std::map<key_type, int> tmp;

        for (const auto& item : recursed) {
            //        for (auto it=recursed.begin(); it!=recursed.end(); ++it) {
            tmp[ftmul(k1l, item.first)[0].first] += item.second;
            //            tmp[ftmul(k1l, it->first)[0].first] += it->second;
        }
        return {tmp.begin(), tmp.end()};
    }

public:
    result_type operator()(argument_type lhs, argument_type rhs) const
    {
        static const boost::container::small_vector<pair_type, 0> null;

        auto lhs_deg = lhs.size();
        if (lhs_deg == 0) {
            return null;
        }
        if (Depth > 0 && ((lhs_deg + rhs.size()) > Depth)) {
            return null;
        }

        return base::cached_compute(lhs, rhs);
    }
};

template<DEG Width, DEG Depth>
using half_shuffle_multiplier = left_half_shuffle_multiplier<Width, Depth>;

template<DEG Width, DEG Depth>
class shuffle_tensor_multiplier
    : multiplier_base<shuffle_tensor_multiplier<Width, Depth>, tensor_basis<Width, Depth>>,
      left_half_shuffle_multiplier<Width, Depth>
{
    using base = multiplier_base<shuffle_tensor_multiplier<Width, Depth>, tensor_basis<Width, Depth>>;
    friend base;

    using half_shuffle_base = left_half_shuffle_multiplier<Width, Depth>;

public:
    using basis_type = tensor_basis<Width, Depth>;
    using key_type = typename basis_type::KEY;
    using scalar_type = int;
    using typename base::pair_type;

    using typename base::inner_result_type;
    using typename base::result_type;
    using argument_type = const key_type&;

private:
    inner_result_type prod_impl(argument_type lhs, argument_type rhs) const
    {
        if (lhs.size() == 0 && rhs.size() == 0) {
            return {{key_type(lhs), 1}};
        }
        return half_shuffle_base::shuffle(lhs, rhs);
    }

public:
#pragma clang diagnostic push
#pragma ide diagnostic ignored "HidingNonVirtualFunction"
    result_type operator()(argument_type lhs, argument_type rhs) const
    {
        static const boost::container::small_vector<pair_type, 0> null;

        if ((lhs.size() + rhs.size()) > Depth) {
            return null;
        }

        return base::cached_compute(lhs, rhs);
        //        return half_shuffle_base::shuffle(lhs, rhs);
    }
#pragma clang diagnostic pop
};

template<DEG Width, DEG Depth>
class area_tensor_multiplier : multiplier_base<area_tensor_multiplier<Width, Depth>, tensor_basis<Width, Depth>>,
                               left_half_shuffle_multiplier<Width, Depth>
{
    using base = multiplier_base<area_tensor_multiplier, tensor_basis<Width, Depth>>;
    using half_shuffle_base = left_half_shuffle_multiplier<Width, Depth>;
    friend base;

public:
    using basis_type = tensor_basis<Width, Depth>;
    using key_type = typename basis_type::KEY;
    using typename base::argument_type;
    using typename base::inner_result_type;
    using typename base::pair_type;
    using typename base::result_type;

private:
    inner_result_type prod_impl(argument_type lhs, argument_type rhs) const
    {
        assert(lhs.valid() && rhs.valid());
        return base::sub(
                half_shuffle_base::operator()(rhs, lhs),
                half_shuffle_base::operator()(lhs, rhs));
    }

public:
#pragma clang diagnostic push
#pragma ide diagnostic ignored "HidingNonVirtualFunction"
    result_type operator()(argument_type lhs, argument_type rhs) const
    {
        static const boost::container::small_vector<pair_type, 0> null;

        if (Depth > 0 && (lhs.size() + rhs.size()) > Depth) {
            return null;
        }

        return base::cached_compute(lhs, rhs);
    }
#pragma clang diagnostic pop
};

template<DEG Width, DEG Depth>
using left_half_shuffle_multiplication =
        base_multiplication<left_half_shuffle_multiplier<Width, Depth>>;

template<DEG Width, DEG Depth>
using half_shuffle_multiplication =
        left_half_shuffle_multiplication<Width, Depth>;

template<DEG Width, DEG Depth>
using area_tensor_multiplication =
        base_multiplication<area_tensor_multiplier<Width, Depth>>;

template<DEG Width, DEG Depth>
using shuffle_tensor_multiplication =
        base_multiplication<shuffle_tensor_multiplier<Width, Depth>>;

template<typename Coeff, DEG n_letters, DEG max_degree, typename...>
class shuffle_tensor;

namespace dtl {

template<typename Coeff, DEG n_letters, DEG max_degree,
         template<typename, typename, typename...> class VectorType,
         typename Derived,
         typename... Args>
class free_tensor_base : public algebra<
                                 free_tensor_basis<n_letters, max_degree>,
                                 Coeff,
                                 free_tensor_multiplication<n_letters, max_degree>,
                                 VectorType,
                                 Derived,
                                 Args...>
{
    using algebra_type = algebra<
            free_tensor_basis<n_letters, max_degree>,
            Coeff,
            free_tensor_multiplication<n_letters, max_degree>,
            VectorType,
            Derived,
            Args...>;

public:
    using basis_type = free_tensor_basis<n_letters, max_degree>;
    using key_type = typename basis_type::KEY;
    using scalar_type = typename Coeff::S;

    // Inherit constructors
    using algebra_type::algebra_type;

    explicit free_tensor_base(typename boost::call_traits<scalar_type>::param_type s)
        : algebra_type(key_type{}, s)
    {}

    template<typename Letter, typename Scalar>
    explicit free_tensor_base(Letter let, Scalar sca)
        : algebra_type(algebra_type::basis.keyofletter(LET(let)), scalar_type(sca))
    {}
};

}// namespace dtl

/**
 * @brief A specialisation of the algebra class with a free tensor basis.
 *
 * Mathematically, the algebra of free_tensor instances is a free associative
 * algebra. With respect to the inherited algebra class, the essential
 * distinguishing feature of this class is the basis class used, and in
 * particular the basis::prod() member function. Thus, the most important
 * information is in the definition of free_tensor_basis. Notice that this
 * associative algebra of free tensors includes as a sub-algebra the
 * associative algebra corresponding to the SCALAR type. This is permitted by
 * the existence of empty keys in free_tensor_basis.
 */
template<typename Coeff, DEG n_letters, DEG max_degree,
         template<typename, typename, typename...> class VectorType,
         typename... Args>
class free_tensor
    : public dtl::free_tensor_base<Coeff,
                                   n_letters,
                                   max_degree,
                                   VectorType,
                                   free_tensor<Coeff, n_letters, max_degree, VectorType, Args...>,
                                   Args...>
{
    typedef free_tensor_multiplication<n_letters, max_degree> multiplication_t;

    using base = dtl::free_tensor_base<Coeff, n_letters, max_degree, VectorType,
                                       free_tensor<Coeff, n_letters, max_degree, VectorType, Args...>, Args...>;

public:
    /// The basis type.
    typedef free_tensor_basis<n_letters, max_degree> BASIS;
    /// Import of the KEY type.
    typedef typename BASIS::KEY KEY;
    /// The algebra type.
    typedef algebra<BASIS, Coeff, multiplication_t, VectorType, free_tensor<Coeff, n_letters, max_degree, VectorType, Args...>, Args...> ALG;

    typedef typename Coeff::SCA SCA;
    typedef typename Coeff::RAT RAT;

    /// The sparse_vector type.
    typedef typename ALG::VECT VECT;

    /// Import of the iterator type.
    typedef typename ALG::iterator iterator;
    /// Import of the constant iterator type.
    typedef typename ALG::const_iterator const_iterator;

public:
    using base::base;

    free_tensor& operator=(const free_tensor&) = default;

    friend free_tensor exp(const free_tensor& arg)
    {
        // Computes the truncated exponential of arg
        // 1 + arg + arg^2/2! + ... + arg^n/n! where n = max_degree
        typename tensor_basis<n_letters, max_degree>::KEY kunit;
        free_tensor result(kunit);
        free_tensor unit(kunit);
        for (DEG i = max_degree; i >= 1; --i) {
            result.mul_scal_div(arg, typename Coeff::Q(i));
            result += unit;
        }
        return result;
    }

    /**
     * Fused multiply exponential operation for free tensors.
     *
     * Computes a*exp(x) using a modified Horner's method for the case when x does
     * not have a constant term. If the argument exp_arg has a constant term, it
     * is ignored.
     *
     * For a real number x, one can expand exp(x) up to degree n as
     *
     *     1 + b_1 x(1 + b_2 x(1 + ... b_n x(1)) ...)
     *
     * where each b_i has the value 1/i. This formulae works when x is a free
     * tensor, or indeed any element in an unital (associative) algebra. Working
     * through the result of multiplying on the left by another element a in the
     * above gives the expansion
     *
     *     a + b1 (a + b_2 (a + ... b_n (a)x) ... x)x.
     *
     * This is the result of a*exp(x). In a non-commutative algebra this need not
     * be equal to exp(x)*a.
     *
     * @param exp_arg free_tensor (const reference) to expentiate (x).
     * @return free_tensor containing a*exp(x)
     */
    free_tensor fmexp(const free_tensor& exp_arg) const
    {
        free_tensor result(*this), x(exp_arg);
        typename free_tensor::KEY kunit;

        auto unit_elt = x.find(kunit);
        if (unit_elt != x.end() && unit_elt->value() != typename free_tensor::SCALAR(0)) {
            x.erase(unit_elt);
        }

        for (DEG i = max_degree; i >= 1; --i) {
            result.mul_scal_div(x, typename free_tensor::SCALAR(i), max_degree - i + 1);
            result += *this;
        }

        return result;
    }

    /// Inplace version of fmexp
    free_tensor& fmexp_inplace(const free_tensor& exp_arg)
    {
        free_tensor original(*this), x(exp_arg);
        typename free_tensor::KEY kunit;
        auto unit_elt = x.find(kunit);

        if (unit_elt != x.end() && unit_elt->value() != typename free_tensor::SCALAR(0)) {
            x.erase(unit_elt);
        }

        for (DEG i = max_degree; i >= 1; --i) {
            this->mul_scal_div(x, typename free_tensor::SCALAR(i), max_degree - i + 1);
            *this += original;
        }

        return *this;
    }

private:
    // Implementation of the antipode for sparse vector types.
    free_tensor antipode_impl(vectors::dtl::access_type_sparse) const
    {
        free_tensor result;

        for (auto cit = this->begin(); cit != this->end(); ++cit) {
            KEY temp_key = cit->key();
            KEY temp_key_reverse = temp_key.reverse();
            SCA temp_value = cit->value();

            int sign;

            if (temp_key.size() % 2 == 0) {
                sign = 1;
            }
            else {
                sign = -1;
            }

            result.add_scal_prod(temp_key_reverse, sign * temp_value);
        }

        return result;
    }

    // Implementation of the antipode for dense vector types.
    free_tensor antipode_impl(vectors::dtl::access_type_dense) const
    {
        free_tensor result;

#ifdef LIBALGEBRA_MAX_TILE_LETTERS
        constexpr DEG CalcLetters = integer_maths::logN(static_cast<unsigned>(LIBALGEBRA_L1_CACHE_SIZE), n_letters) / 2;
        constexpr DEG BlockLetters = (CalcLetters > LIBALGEBRA_MAX_TILE_LETTERS) ? LIBALGEBRA_MAX_TILE_LETTERS : CalcLetters;
#else
        constexpr DEG BlockLetters = integer_maths::logN(LIBALGEBRA_L1_CACHE_SIZE / sizeof(SCA), n_letters) / 2;
#endif

        const auto curr_degree = this->degree();

        vectors::dtl::vector_base_access::convert(result).resize_to_degree(curr_degree);

        // Get the pointers to the start of the data blob in memory.
        const SCA* src_ptr = vectors::dtl::data_access<VectorType<BASIS, Coeff>>::range_begin(vectors::dtl::vector_base_access::convert(*this));
        SCA* dst_ptr = vectors::dtl::data_access<VectorType<BASIS, Coeff>>::range_begin(vectors::dtl::vector_base_access::convert(result));

        dtl::tiled_inverse_operator<n_letters, max_degree, BlockLetters, SCA, dtl::default_signer> t;

        t(src_ptr, dst_ptr, curr_degree);

        return result;
    }

public:
    /// Computes the antipode of a free_tensor instance.
    inline friend free_tensor antipode(const free_tensor& arg)
    {
        // Get the trait to access the storage tag, although it now occurs to me that we already had
        // the vector type as a template argument, so we might be able to dispatch off that instead
        // of a tag. But the tag will do for now.
        using trait = vectors::dtl::data_access<VectorType<BASIS, Coeff>>;

        // Now use tagged dispatch to pick the correct implementation
        return arg.antipode_impl(typename trait::tag());
    }

#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
private:
    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive& ar, unsigned int const /* version */)
    {
        ar& boost::serialization::base_object<base>(*this);
    }
#endif
};

/*

/// Computes the truncated exponential of a free_tensor instance.
template<typename Coeffs, DEG Width, DEG Depth, template<typename, typename, typename...> class VectorType, typename... Args>
free_tensor<Coeffs, Width, Depth, VectorType, Args...>
exp(const free_tensor<Coeffs, Width, Depth, VectorType, Args...>& arg)
{
    // Computes the truncated exponential of arg
    // 1 + arg + arg^2/2! + ... + arg^n/n! where n = max_degree
    typename tensor_basis<Width, Depth>::KEY kunit;
    free_tensor<Coeffs, Width, Depth, VectorType, Args...> result(kunit);
    free_tensor<Coeffs, Width, Depth, VectorType, Args...> unit(kunit);
    for (DEG i = Depth; i >= 1; --i) {
        result.mul_scal_div(arg, typename Coeffs::Q(i));
        result += unit;
    }
    return result;
}
*/

/// Computes the truncated logarithm of a free_tensor instance.
template<typename Coeffs, DEG Width, DEG Depth, template<typename, typename, typename...> class VectorType, typename... Args>
free_tensor<Coeffs, Width, Depth, VectorType, Args...>
log(const free_tensor<Coeffs, Width, Depth, VectorType, Args...>& arg)
{
    // Computes the truncated log of arg up to degree max_degree
    // The coef. of the constant term (empty word in the monoid) of arg
    // is forced to 1.
    // log(arg) = log(1+x) = x - x^2/2 + ... + (-1)^(n+1) x^n/n.
    // max_degree must be > 0
    using ft_type = free_tensor<Coeffs, Width, Depth, VectorType, Args...>;

    typename tensor_basis<Width, Depth>::KEY kunit;
    ft_type tunit(kunit);
    ft_type x(arg);
    auto it = x.find(kunit);
    if (it != x.end()) {
        x.erase(it);
    }
    ft_type result;

    for (DEG i = Depth; i >= 1; --i) {
        if (i % 2 == 0) {
            result.sub_scal_div(tunit, typename Coeffs::Q(i));
        }
        else {
            result.add_scal_div(tunit, typename Coeffs::Q(i));
        }
        result *= x;
    }

    return result;
}

/// Computes the truncated inverse of a free_tensor instance.
template<typename Coeffs, DEG Width, DEG Depth, template<typename, typename, typename...> class VectorType, typename... Args>
free_tensor<Coeffs, Width, Depth, VectorType, Args...>
inverse(const free_tensor<Coeffs, Width, Depth, VectorType, Args...>& arg)
{
    // Computes the truncated inverse of arg up to degree max_degree
    // An exception is thrown if the leading term is zero.
    // the module assumes
    // (a+x)^(-1) = (a(1+x/a))^(-1)
    //  = a^(-1)(1 - x/a + x^2/a^2 + ... + (-1)^(n) x^n/a^n)
    // = a^(-1) - x/a*[a^(-1)(1 - x/a + x^2/a^2 + ... + (-1)^(n)
    // x^(n-1)/a^(n-1)))]. S_n = a^(-1) + z S_{n-1}; z = - x/a ; S_0 = a^(-1)
    // max_degree must be > 0
    using ft_type = free_tensor<Coeffs, Width, Depth, VectorType, Args...>;

    typename tensor_basis<Width, Depth>::KEY kunit;
    typename Coeffs::S a(0);
    ft_type x, z(a);

    auto it(arg.find(kunit));
    if (it == arg.end()) {
        // const term a is 0;
        throw std::invalid_argument("divide-by-zero");
    }
    else {
        a = (*it).value();
        x = arg;
        x.erase(kunit);
    }

    // S_n = a + z S_{ n - 1 }; z = -x / a; S_0 = a
    //
    //  the nonzero scalar component a of the tensor arg restored to a tensor
    ft_type free_tensor_a_inverse(typename Coeffs::S(1) / a), result(free_tensor_a_inverse);
    // z := - x/a
    z.sub_scal_div(x, a);
    // the iteration
    for (DEG i = 0; i != Depth; ++i) {
        auto tmp = z*result;
        result = free_tensor_a_inverse + z * result;
    }
    return result;
}

/// Computes the truncated inverse of a free_tensor instance.
template<typename Coeffs, DEG Width, DEG Depth, template<typename, typename, typename...> class VectorType, typename... Args>
free_tensor<Coeffs, Width, Depth, VectorType, Args...>
reflect(const free_tensor<Coeffs, Width, Depth, VectorType, Args...>& arg)
{
    return antipode(arg);
}

/**
 * @brief A specialisation of the algebra class with a shuffle tensor basis.
 *
 * Mathematically, the algebra of shuffle_tensor instances is a shuffle
 * associative algebra associated to a free associative algebra. With respect
 * to the inherited algebra class, the essential distinguishing feature of
 * this class is the basis class used, and in particular the basis::prod()
 * member function. Thus, the most important information is in the definition
 * of shuffle_tensor_basis. Notice that this associative algebra of free
 * tensors includes as a sub-algebra the associative algebra corresponding to
 * the SCALAR type. This is permitted by the existence of empty keys in
 * shuffle_tensor_basis.
 */
template<typename Coeff, DEG n_letters, DEG max_degree, typename...>
class shuffle_tensor : public algebra<
                               shuffle_tensor_basis<n_letters, max_degree>,
                               Coeff,
                               shuffle_tensor_multiplication<n_letters, max_degree>>
{
    typedef shuffle_tensor_multiplication<n_letters, max_degree> multiplication_t;

public:
    /// The basis type.
    typedef shuffle_tensor_basis<n_letters, max_degree> BASIS;
    /// Import of the KEY type.
    typedef typename BASIS::KEY KEY;
    /// The algebra type.
    typedef algebra<BASIS, Coeff, multiplication_t> ALG;

    typedef typename Coeff::SCA SCA;
    typedef typename Coeff::RAT RAT;

    /// The sparse_vector type.
    typedef typename ALG::VECT VECT;

    /// Import of the iterator type.
    typedef typename ALG::iterator iterator;
    /// Import of the constant iterator type.
    typedef typename ALG::const_iterator const_iterator;

public:
    using ALG::ALG;

    /// Default constructor.
    shuffle_tensor()
    {}

    /// Copy constructor.
    shuffle_tensor(const shuffle_tensor& t)
        : ALG(t)
    {}

    /// Constructs an instance from a free_tensor instance.
    template<template<typename, typename, typename...> class VectorType, typename... Args>
    shuffle_tensor(const free_tensor<Coeff, n_letters, max_degree, VectorType, Args...>& t)
    {
        typename free_tensor<Coeff, n_letters, max_degree, VectorType, Args...>::const_iterator i;
        for (i = t.begin(); i != t.end(); ++i) {
            (*this)[i->key()] += i->value();
        }
    }

    /// Constructs an instance from an algebra instance.
    shuffle_tensor(const ALG& a)
        : ALG(a)
    {}

    /// Constructs an instance from a sparse_vector instance.
    shuffle_tensor(const VECT& v)
        : ALG(v)
    {}

    /// Constructs a unidimensional instance from a letter and a scalar.
    shuffle_tensor(LET letter, const SCA& s)
        : ALG(VECT::basis.keyofletter(letter), s)
    {}

    /// Constructs a unidimensional instance from a key (basis element).
    explicit shuffle_tensor(const KEY& k)
        : ALG(k)
    {}

    /// Constructs a unidimensional instance from a scalar.
    explicit shuffle_tensor(const SCA& s)
        : ALG(VECT::basis.empty_key, s)
    {}

    shuffle_tensor& operator=(const shuffle_tensor&) = default;
    shuffle_tensor& operator=(shuffle_tensor&&) noexcept = default;

#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
private:
    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive& ar, unsigned int const /* version */)
    {
        ar& boost::serialization::base_object<ALG>(*this);
    }
#endif
};

template<typename Coeff, DEG Width, DEG Depth,
         template<typename, typename, typename...> class VectorType,
         typename... Args>
using half_shuffle_tensor = algebra<half_shuffle_tensor_basis<Width, Depth>,
                                    Coeff,
                                    half_shuffle_multiplication<Width, Depth>,
                                    VectorType,
                                    Args...>;

template<typename Coeff, DEG Width, DEG Depth,
         template<typename, typename, typename...> class VectorType,
         typename... Args>
using area_tensor = algebra<area_tensor_basis<Width, Depth>,
                            Coeff,
                            area_tensor_multiplication<Width, Depth>,
                            VectorType,
                            Args...>;

namespace vectors {

namespace dtl {

#define LIBALGEBRA_DENSE_TENSOR_REF_BINOP(OP)         \
    template<typename T>                              \
    dense_tensor_item_reference& operator OP(T other) \
    {                                                 \
        *m_item OP scalar_type(other);                \
        if (m_reverse_item != nullptr) {              \
            *m_reverse_item OP scalar_type(other);    \
        }                                             \
        return *this;                                 \
    }

template<typename Coeffs>
class dense_tensor_item_reference
{
    using scalar_type = typename Coeffs::S;
    scalar_type* m_item;
    scalar_type* m_reverse_item;

public:
    explicit dense_tensor_item_reference(scalar_type* item, scalar_type* reverse_item) noexcept
        : m_item(item), m_reverse_item(reverse_item)
    {
        assert(m_reverse_item == nullptr || *m_item == *m_reverse_item);
    }

    operator const scalar_type&() noexcept { return *m_item; }

    template<typename T>
    dense_tensor_item_reference& operator=(const T& new_val)
    {
        *m_item = scalar_type(new_val);
        if (m_reverse_item != nullptr) {
            *m_reverse_item = scalar_type(new_val);
        }
        return *this;
    }

    template<typename T>
    dense_tensor_item_reference& operator=(T&& new_val) noexcept
    {
        *m_item = scalar_type(std::forward<T>(new_val));
        if (m_reverse_item != nullptr) {
            *m_reverse_item = *m_item;
        }
        return *this;
    }

    LIBALGEBRA_DENSE_TENSOR_REF_BINOP(+=)
    LIBALGEBRA_DENSE_TENSOR_REF_BINOP(-=)
    LIBALGEBRA_DENSE_TENSOR_REF_BINOP(*=)
    LIBALGEBRA_DENSE_TENSOR_REF_BINOP(/=)
    LIBALGEBRA_DENSE_TENSOR_REF_BINOP(&=)
    LIBALGEBRA_DENSE_TENSOR_REF_BINOP(|=)
};

#undef LIBALGEBRA_DENSE_TENSOR_REF_BINOP

template<DEG Width, DEG Depth, typename Coeffs>
class dense_tensor_const_iterator_item
{
    using scalar_type = typename Coeffs::S;
    const scalar_type* m_ptr;

    using vector_type = dense_vector<free_tensor_basis<Width, Depth>, Coeffs>;

    friend class iterators::vector_iterator<dense_tensor_const_iterator_item>;
    friend vector_type;

    const vector_type* p_vector;

public:
    using basis_type = free_tensor_basis<Width, Depth>;
    using key_type = typename free_tensor_basis<Width, Depth>::KEY;
    using reference = const scalar_type&;

    dense_tensor_const_iterator_item() = default;

    explicit dense_tensor_const_iterator_item(const vector_type& vect,
                                              const scalar_type* item)
        : m_ptr(item), p_vector(&vect)
    {}

    DIMN index() const
    {
        return static_cast<DIMN>(m_ptr - p_vector->as_ptr());
    }

    key_type key()
    {
        assert(index() < p_vector->dimension());
        return basis_type::index_to_key(index());
    }

    reference value()
    {
        return *m_ptr;
    }

private:
    bool compare_iterators(const dense_tensor_const_iterator_item& other) const
    {
        if (m_ptr == nullptr || other.m_ptr == nullptr) {
            return true;
        }
        return m_ptr >= other.m_ptr;
    }

    void advance()
    {
        if (m_ptr == nullptr) {
            return;
        }
        const auto* end = p_vector->as_ptr() + p_vector->dimension();
        do {
            ++m_ptr;
        } while (m_ptr != end && *m_ptr == Coeffs::zero);
    }
};

template<DEG Width, DEG Depth, typename Coeffs>
class dense_tensor_iterator_item
{
    using tsi = alg::dtl::tensor_size_info<Width>;
    using vector_type = dense_vector<free_tensor_basis<Width, Depth>, Coeffs>;
    using scalar_type = typename Coeffs::S;

    scalar_type* m_ptr;
    vector_type* p_vector;
    DEG m_degree = 0;

    friend class iterators::vector_iterator<dense_tensor_iterator_item>;
    friend vector_type;

public:
    using basis_type = free_tensor_basis<Width, Depth>;
    using key_type = typename free_tensor_basis<Width, Depth>::KEY;
    using reference = scalar_type&;

    dense_tensor_iterator_item() = default;

    explicit dense_tensor_iterator_item(vector_type& vect,
                                        scalar_type* item)
        : m_ptr(item),
          p_vector(&vect)
    {
        auto it = std::lower_bound(tsi::degree_sizes.begin(), tsi::degree_sizes.end(), index(), std::less_equal<>());
        m_degree = static_cast<DEG>(it - tsi::degree_sizes.begin());
    }

    DIMN index() const
    {
        return static_cast<DIMN>(m_ptr - p_vector->as_mut_ptr());
    }

    key_type key()
    {
        assert(index() < p_vector->dimension());
        return basis_type::index_to_key(index());
    }

    reference value()
    {
        return *m_ptr;
    }

private:
    bool compare_iterators(const dense_tensor_iterator_item& other) const
    {
        if (m_ptr == nullptr || other.m_ptr == nullptr) {
            return true;
        }
        return m_ptr >= other.m_ptr;
    }

    void advance()
    {
        auto* end = p_vector->as_mut_ptr() + p_vector->dimension();
        do {
            ++m_ptr;
        } while (m_ptr != end && *m_ptr == Coeffs::zero);
    }
};

}// namespace dtl
template<DEG Width, DEG Depth, typename Coeffs>
class dense_vector<free_tensor_basis<Width, Depth>, Coeffs>
    : protected base_vector<free_tensor_basis<Width, Depth>, Coeffs>
{
    using storage_type = dense_storage<typename Coeffs::S>;
    using base_vector_type = base_vector<free_tensor_basis<Width, Depth>, Coeffs>;
    using tsi = ::alg::dtl::tensor_size_info<Width>;
    friend class dtl::data_access_base<dense_vector>;
    friend class dtl::dense_tensor_iterator_item<Width, Depth, Coeffs>;

    friend class traditional_free_tensor_multiplication<Width, Depth>;


    class reverse_storage : public storage_type
    {
        bool m_valid = false;

    public:

        constexpr bool valid() const noexcept
        { return m_valid; }

        constexpr operator bool() const noexcept
        { return valid() && storage_type::size() > 0; }

        constexpr void validate() noexcept { m_valid = true; }
        constexpr void invalidate() noexcept { m_valid = false; }


    };


    storage_type m_data;
//    std::vector<typename Coeffs::S> m_reverse_data;
    reverse_storage m_reverse_data;
    DIMN m_dimension = 0;
    DEG m_degree = 0;

    using basis_traits = basis::basis_traits<free_tensor_basis<Width, Depth>>;

    static typename tsi::degree_sizes_t::const_iterator index_to_size_iter(DIMN index) noexcept
    {
        return std::upper_bound(tsi::degree_sizes.begin(), tsi::degree_sizes.end(), index, std::less_equal<>());
    }

    static DEG size_iter_to_deg(typename tsi::degree_sizes_t::const_iterator it) noexcept
    {
        return static_cast<DEG>(it - tsi::degree_sizes.begin());
    }

    static DIMN reverse_index(DIMN index_in_degree, DEG deg) noexcept
    {
        DIMN result = 0;
        for (DEG i = 0; i < deg; ++i) {
            result *= DIMN(Width);
            result += index_in_degree % Width;
            index_in_degree /= Width;
        }

        return result;
    }

    void construct_reverse_data(DEG degree)
    {
#ifdef LIBALGEBRA_MAX_TILE_LETTERS
        constexpr DEG CalcLetters = integer_maths::logN(static_cast<unsigned>(LIBALGEBRA_L1_CACHE_SIZE), n_letters) / 2;
        constexpr DEG BlockLetters = (CalcLetters > LIBALGEBRA_MAX_TILE_LETTERS) ? LIBALGEBRA_MAX_TILE_LETTERS : CalcLetters;
#else
        constexpr DEG BlockLetters = ::alg::integer_maths::logN(LIBALGEBRA_L1_CACHE_SIZE / sizeof(scalar_type), Width) / 2;
#endif
        const auto& dsizes = tsi::degree_sizes;

        if (degree > 0) {
            if (m_reverse_data.size() < tsi::degree_sizes[degree]) {
                m_reverse_data.resize(tsi::degree_sizes[degree]);
            }
            assert(m_reverse_data.size() >= tsi::degree_sizes[degree]);
            ::alg::dtl::tiled_inverse_operator<Width, Depth, BlockLetters, scalar_type, ::alg::dtl::non_signing_signer> t;
            t(m_data.begin(), m_reverse_data.data(), degree);
            assert(*m_data.begin() == *m_reverse_data.begin());
        }
        m_reverse_data.validate();
    }

public:
    using basis_type = free_tensor_basis<Width, Depth>;
    using key_type = typename basis_type::KEY;
    using coefficient_ring = Coeffs;
    using scalar_type = typename Coeffs::S;
    using rational_type = typename Coeffs::Q;

    using scalar_param_type = typename boost::call_traits<scalar_type>::param_type;

    using pointer = typename storage_type::pointer;
    using const_pointer = typename storage_type::const_pointer;
    using reference = scalar_type&;
    using const_reference = typename storage_type::const_reference;

    using iterator = iterators::vector_iterator<
            dtl::dense_tensor_iterator_item<Width, Depth, Coeffs>>;
    using const_iterator = iterators::vector_iterator<
            dtl::dense_tensor_const_iterator_item<Width, Depth, Coeffs>>;

    using base_vector_type::degree_tag;

    using BASIS = basis_type;
    using COEFFS = coefficient_ring;
    using KEY = key_type;
    using SCALAR = scalar_type;
    using RATIONAL = rational_type;

    using base_vector_type::basis;
    using base_vector_type::mone;
    using base_vector_type::one;
    using base_vector_type::zero;

    dense_vector() = default;
    dense_vector(const dense_vector&) = default;
    dense_vector(dense_vector&&) noexcept = default;

    explicit dense_vector(key_type key, scalar_param_type scalar = Coeffs::one)
        : m_data(), m_reverse_data(), m_dimension(0), m_degree(0)
    {
        auto idx = basis_traits::key_to_index(base_vector_type::basis, key);
        resize_to_dimension(idx+1);
        m_data[idx] = scalar;
    }

    dense_vector(const_pointer begin, const_pointer end)
        : m_data(), m_reverse_data(), m_dimension(0), m_degree(0)
    {
        resize_to_dimension(static_cast<DIMN>(end - begin));
        std::copy(begin, end, m_data.begin());
//        if (m_degree > 0) {
//            construct_reverse_data(m_degree-1);
//        }
    }

    dense_vector& operator=(const dense_vector&) = default;
    dense_vector& operator=(dense_vector&&) noexcept = default;

    const_reference value(DIMN index) const noexcept
    {
        return m_data[index];
    }

    reference value(DIMN index) noexcept
    {
//        resize_to_dimension(index + 1);
        m_reverse_data.invalidate();
        return m_data[index];
    }

    const_reference operator[](key_type key) const noexcept
    {
        auto idx = basis_traits::key_to_index(base_vector_type::basis, key);
        if (idx < m_dimension) {
            return m_data[idx];
        }
        return Coeffs::zero;
    }

    reference operator[](key_type key) noexcept
    {
        auto idx = resize_for_key(key);
        assert(key.size() <= m_degree);
        return value(idx);
    }

    /// Reserve to dimension
    void reserve_to_dimension(const DIMN dim)
    {
        if (dim > m_dimension) {
            auto info = basis_traits::next_resize_dimension(base_vector_type::basis, dim, m_degree);
            m_data.reserve(info.size);
            m_dimension = info.dimension;
            m_degree = info.degree;
            m_reverse_data.invalidate();
        }
        assert(m_data.size() == m_dimension);
    }

    /// Reserve to degree
    void reserve_to_degree(const DEG deg)
    {
        if (deg > m_degree) {
            auto info = basis_traits::next_resize_dimension(base_vector_type::basis, 0, deg);
            m_data.reserve(info.size);
            m_dimension = info.dimension;
            m_degree = info.degree;
            m_reverse_data.invalidate();
        }
        assert(m_data.size() == m_dimension);

    }

    DIMN resize_for_key(key_type key)
    {
        auto info = basis_traits::key_resize_dimension(base_vector_type::basis, key);
        if (info.dimension > m_dimension) {
            m_data.resize(info.size);
            m_dimension = info.dimension;
            m_degree = info.degree;
            m_reverse_data.invalidate();
        }
        assert(m_data.size() == m_dimension);
        return basis_traits::key_to_index(base_vector_type::basis, key);
    }

    /// Resize to dimension
    void resize_to_dimension(DIMN dim)
    {
        auto info = basis_traits::next_resize_dimension(base_vector_type::basis, dim, m_degree);
        if (info.dimension > m_dimension) {
            m_data.resize(info.size);
            m_dimension = info.dimension;
            m_degree = info.degree;
            m_reverse_data.invalidate();
        }
        else if (info.dimension == m_dimension) {
            m_degree = info.degree;
        }
        assert(m_data.size() == m_dimension);
    }

    /// Reserve to degree
    void resize_to_degree(const DEG deg)
    {
        if (deg > m_degree || m_dimension == 0) {
            auto info = basis_traits::next_resize_dimension(base_vector_type::basis, 0, deg);
            m_data.resize(info.size);
            m_dimension = info.dimension;
            m_degree = info.degree;
            m_reverse_data.invalidate();
        }
    }

    DEG degree() const noexcept { return m_degree; }
    DIMN dimension() const noexcept { return m_dimension; };
    DIMN reverse_dimension() const noexcept { return m_reverse_data.size(); }
    bool degree_equals(DEG degree) const noexcept { return m_degree == degree; }

    pointer as_mut_ptr() noexcept { return m_data.begin(); }
    const_pointer as_ptr() const noexcept { return m_data.begin(); }
    const_pointer as_rptr() const noexcept { return m_reverse_data.data(); }
    pointer as_mut_rptr() noexcept { return m_reverse_data.data(); }

    reverse_storage& reverse_data() noexcept { return m_reverse_data; }

protected:
    /// Index of key and key of index
    static KEY index_to_key(const DIMN index)
    {
        return basis.index_to_key(index);
    }

    /// Get the index corresponding to a key
    static DIMN key_to_index(const KEY& key)
    {
        return basis.key_to_index(key);
    }

    /// Get the degree of the key at an index
    static DEG index_to_degree(const DIMN idx)
    {
        return basis.degree(index_to_key(idx));
    }

    /// Get the index at which the elements of degree start
    static DIMN start_of_degree(const DEG deg)
    {
        return basis.start_of_degree(deg);
    }

    /// Get he maximum feasible dimension for bases with degree
    template<DEG D>
    static DIMN max_dimension(alg::basis::with_degree<D>)
    {
        return basis.start_of_degree(degree_tag.max_degree + 1);
    }

    /// Get the maximum feasible dimension for a basis without degree
    static DIMN max_dimension(alg::basis::without_degree)
    {
        return basis.max_dimension();
    }

public:
    bool empty() const noexcept
    {
        if (m_data.empty()) {
            return true;
        }

        for (DIMN i = 0; i < m_dimension; ++i) {
            if (m_data[i] != Coeffs::zero) {
                return false;
            }
        }
        return true;
    }

    DIMN size() const noexcept
    {
        DIMN result = 0;
        for (DIMN i = 0; i < m_dimension; ++i) {
            result += static_cast<DIMN>(m_data[i] != Coeffs::zero);
        }
        return result;
    }

    void swap(dense_vector& other)
    {
        m_data.swap(other.m_data);
        m_reverse_data.swap(other.m_reverse_data);
        std::swap(m_dimension, other.m_dimension);
        std::swap(m_degree, other.m_degree);
    }

    template<typename Scalar>
    void insert(key_type key, Scalar val)
    {
        operator[](key) = val;
    }

    void erase(key_type key)
    {
        operator[](key) = Coeffs::zero;
    }

    void erase(iterator& it)
    {
        it->value() = Coeffs::zero;

        if (it->index() < m_reverse_data.size() && m_reverse_data) {
            m_reverse_data[basis_type::key_to_index(it->key().reverse())] = Coeffs::zero;
        }

    }

    void clear()
    {
        std::fill(m_data.begin(), m_data.end(), Coeffs::zero);
        std::fill(m_reverse_data.begin(), m_reverse_data.end(), Coeffs::zero);
    }

    iterator begin()
    {
        m_reverse_data.invalidate();
        return iterator(*this, m_data.begin());
    }
    iterator end()
    {
        m_reverse_data.invalidate();
        return iterator(*this, m_data.end());
    }
    const_iterator begin() const
    {
        return const_iterator(*this, m_data.begin());
    }
    const_iterator end() const
    {
        return const_iterator(*this, m_data.end());
    }
    const_iterator cbegin() const
    {
        return const_iterator(*this, m_data.begin());
    }
    const_iterator cend() const
    {
        return const_iterator(*this, m_data.end());
    }

    const_iterator find(key_type key) const
    {
        auto deg = key.size();
        if (deg <= m_degree) {
            return const_iterator(*this, m_data.begin() + basis_type::key_to_index(key));
        }
        return end();
    }
    iterator find(key_type key)
    {
        m_reverse_data.invalidate();
        auto deg = key.size();
        if (deg <= m_degree) {
            return iterator(*this, m_data.begin() + basis_type::key_to_index(key));
        }
        return end();
    }
    const_iterator find_index(DIMN idx) const
    {
        if (idx < m_dimension) {
            return const_iterator(*this, m_data.begin() + idx);
        }
        return end();
    }
    iterator find_index(DIMN idx)
    {
        m_reverse_data.invalidate();
        if (idx < m_dimension) {
            return iterator(*this, m_data.begin() + idx);
        }
        return end();
    }

    template<typename InputIt>
    void insert(InputIt begin, InputIt end)
    {
        m_reverse_data.invalidate();
        for (auto it = begin; it != end; ++it) {
            add_scal_prod(it->first, it->second);
        }
//        m_reverse_data.clear();
    }

    template<typename Scalar>
    void add_scal_prod(key_type key, Scalar s)
    {
        m_reverse_data.invalidate();
        operator[](key) += scalar_type(s);
//        m_reverse_data.clear();
    }

    template<typename Fn>
    static void apply_unary_operation(
            dense_vector& result,
            const dense_vector& arg,
            Fn fn)
    {
        assert(result.m_dimension == 0);
        result.reserve_to_dimension(arg.m_dimension);
        for (DIMN i = 0; i < arg.m_dimension; ++i) {
            result.m_data.emplace(i, fn(arg.m_data[i]));
        }
        if (arg.m_reverse_data) {
            result.m_reverse_data.reserve(arg.m_reverse_data.size());
            for (DIMN i=0; i<arg.m_reverse_data.size(); ++i) {
                result.m_reverse_data.emplace(i, fn(arg.m_reverse_data[i]));
            }
        }

    }

    template<typename Fn>
    static void apply_inplace_unary_op(dense_vector& arg, Fn op)
    {
        for (DIMN i = 0; i < arg.m_dimension; ++i) {
            arg.m_data[i] = op(arg.m_data[i]);
        }
        for (DIMN i=0; i<arg.m_reverse_data.size(); ++i) {
            arg.m_reverse_data[i] = op(arg.m_reverse_data[i]);
        }

    }

    template<typename Fn>
    static void apply_flat_binary_operation(
            dense_vector& result,
            const dense_vector& lhs,
            const dense_vector& rhs,
            Fn op)
    {
        auto minmax = std::minmax(lhs.m_dimension, rhs.m_dimension);

        result.reserve_to_dimension(minmax.second);

        const auto size_fwd = minmax.first;

        DIMN i=0;
        for (; i < size_fwd; ++i) {
            result.m_data.emplace(i, op(lhs.m_data[i], rhs.m_data[i]));
        }

        for (; i < lhs.m_dimension; ++i) {
            result.m_data.emplace(i, op(lhs.m_data[i], Coeffs::zero));
        }
        for (; i < rhs.m_dimension; ++i) {
            result.m_data.emplace(i, op(Coeffs::zero, rhs.m_data[i]));
        }
        assert(i == result.m_dimension);

        if (lhs.m_reverse_data && rhs.m_reverse_data) {
            i = 0;
            const auto lhs_rsize = lhs.m_reverse_data.size();
            const auto rhs_rsize = rhs.m_reverse_data.size();

            auto rminmax = std::minmax(lhs_rsize, rhs_rsize);

            result.m_reverse_data.reserve(rminmax.second);

            for (; i<rminmax.first; ++i) {
                result.m_reverse_data.emplace(i, op(lhs.m_reverse_data[i],
                                                    rhs.m_reverse_data[i]));
            }
            for (; i<lhs_rsize; ++i) {
                result.m_reverse_data.emplace(i, op(lhs.m_reverse_data[i],
                                                    Coeffs::zero));
            }
            for (; i<rhs_rsize; ++i) {
                result.m_reverse_data.emplace(i, op(Coeffs::zero,
                                                    rhs.m_reverse_data[i]));
            }
            result.m_reverse_data.validate();
        } else {
            result.m_reverse_data.invalidate();
        }

    }

    template<typename Fn>
    static void apply_inplace_flat_binary_op(
            dense_vector& lhs,
            const dense_vector& rhs,
            Fn op)
    {
        auto old_lhs_dim = lhs.m_dimension;
        auto minmax = std::minmax(old_lhs_dim, rhs.m_dimension);

        const auto lhs_rsize = lhs.m_reverse_data.size();
        const auto rhs_rsize = rhs.m_reverse_data.size();
        const auto size_rev = std::min(lhs_rsize, rhs_rsize);

        if (minmax.second > lhs.m_dimension) {
            lhs.reserve_to_dimension(minmax.second);
        }
        const auto size_fwd = minmax.first;

        DIMN i=0;

        for (; i < size_fwd; ++i) {
            lhs.m_data[i] = op(lhs.m_data[i], rhs.m_data[i]);
        }
        for (; i < old_lhs_dim; ++i) {
            lhs.m_data[i] = op(lhs.m_data[i], Coeffs::zero);
        }
        for (; i < rhs.m_dimension; ++i) {
            lhs.m_data.emplace(i, op(Coeffs::zero, rhs.m_data[i]));
        }
        assert(i == lhs.m_dimension);

        if (lhs.m_reverse_data && rhs.m_reverse_data) {

            const auto lhs_rsize = lhs.m_reverse_data.size();
            const auto rhs_rsize = rhs.m_reverse_data.size();
            const auto rminmax = std::minmax(lhs_rsize, rhs_rsize);

            lhs.m_reverse_data.reserve(rminmax.second);

            i=0;
            for (; i<rminmax.first; ++i) {
                lhs.m_reverse_data[i] = op(lhs.m_reverse_data[i],
                                           rhs.m_reverse_data[i]);
            }
            for (; i<lhs_rsize; ++i) {
                lhs.m_reverse_data[i] = op(lhs.m_reverse_data[i], Coeffs::zero);
            }
            for (; i<rhs_rsize; ++i) {
                lhs.m_reverse_data.emplace(i, op(Coeffs::zero, rhs.m_reverse_data[i]));
            }
            lhs.m_reverse_data.validate();
        } else {
            lhs.m_reverse_data.invalidate();
        }


    }

    scalar_type NormL1() const
    {
        auto ans = Coeffs::zero;
        for (DIMN i = 0; i < m_dimension; ++i) {
            Coeffs::add_inplace(ans, abs(m_data[i]));
        }
        return ans;
    }

    scalar_type NormL1(const DEG degree) const
    {
        auto begin = std::min(m_dimension, basis_type::start_of_degree(degree));
        auto end = std::min(m_dimension, basis_type::start_of_degree(degree + 1));

        auto ans = Coeffs::zero;
        for (DIMN i = begin; i < end; ++i) {
            Coeffs::add_inplace(ans, abs(m_data[i]));
        }
        return ans;
    }

    scalar_type NormLInf() const
    {
        auto ans = Coeffs::zero;
        for (DIMN i = 0; i < m_dimension; ++i) {
            auto abs_val = abs(m_data[i]);
            ans = (abs_val > ans) ? abs_val : ans;
        }
        return ans;
    }

    scalar_type NormLInf(const DEG degree) const
    {
        auto begin = std::min(m_dimension, basis_type::start_of_degree(degree));
        auto end = std::min(m_dimension, basis_type::start_of_degree(degree + 1));

        auto ans = Coeffs::zero;
        for (DIMN i = begin; i < end; ++i) {
            auto abs_val = abs(m_data[i]);
            ans = (abs_val > ans) ? abs_val : ans;
        }
        return ans;
    }

protected:
    std::pair<DIMN, bool> equal_to_min(const dense_vector& rhs) const
    {
        DIMN mid = std::min(m_dimension, rhs.m_dimension);

        for (DIMN i = 0; i < mid; ++i) {
            if (m_data[i] != rhs.m_data[i]) {
                return {mid, false};
            }
        }

        return {mid, true};
    }

    void print_members(std::ostream& os) const
    {
        std::pair<basis_type*, key_type> token;
        token.first = &base_vector_type::basis;
        for (DIMN i = 0; i < m_dimension; ++i) {
            if (m_data[i] != Coeffs::zero) {
                token.second = basis_type::index_to_key(i);
                os << ' ' << m_data[i] << '(' << token << ')';
            }
        }
    }

public:
    bool operator==(const dense_vector& rhs) const
    {
        auto mid_eq = equal_to_min(rhs);
        if (!mid_eq.second) {
            return false;
        }

        for (auto i=mid_eq.first; i<m_dimension; ++i) {
            if (m_data[i] != Coeffs::zero) {
                return false;
            }
        }

        for (auto i=mid_eq.first; i<rhs.m_dimension; ++i) {
            if (rhs.m_data[i] != Coeffs::zero) {
                return false;
            }
        }

        return true;
    }

    bool operator!=(const dense_vector& rhs) const
    {
        return !operator==(rhs);
    }

    friend std::ostream& operator<<(std::ostream& os, const dense_vector& rhs)
    {
        os << '{';
        rhs.print_members(os);
        os << " }";
        return os;
    }

    static bool comp(std::pair<key_type, scalar_type> lhs, std::pair<key_type, scalar_type> rhs)
    {
        return lhs.first < rhs.first;
    }

#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
private:
    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive& ar, const unsigned version)
    {
        ar& m_dimension;
        ar& m_degree;
        ar& m_data;
        ar& m_reverse_data;
    }
#endif
};

}// namespace vectors

}// namespace alg
// Include once wrapper
#endif// DJC_COROPA_LIBALGEBRA_TENSORH_SEEN

// EOF.
