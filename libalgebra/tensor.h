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

#include <boost/container/small_vector.hpp>
#include <boost/functional/hash.hpp>

#include "area_tensor_basis.h"
#include "base_vector.h"
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

    template <typename B>
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
        assert(offset + basis_type::start_of_degree(d)<=basis_type::start_of_degree(d+1));
        return ((lhs_data == nullptr) ? out_data : lhs_data) + offset + basis_type::start_of_degree(d);
    }
    const_pointer right_fwd_read(IDEG d, IDIMN offset = 0) const noexcept
    {
        assert(d >= 0 && d <= rhs_deg);
        assert(offset + basis_type::start_of_degree(d)<=basis_type::start_of_degree(d+1));
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
    const_pointer left_reverse_read_ptr;

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
            reverse_data.resize(basis_type::start_of_degree(base::lhs_deg+1));
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
            reverse_data.resize(basis_type::start_of_degree(base::lhs_deg+1));
            dtl::tiled_inverse_operator<Width, (Depth > 0) ? Depth - 1 : 0, TileLetters, scalar_type, dtl::non_signing_signer> reverser;
            reverser(base::out_data, reverse_data.data(), base::lhs_deg);
            left_reverse_read_ptr = reverse_data.data();
        }
    }

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
        assert(index < static_cast<IDEG>(tsi::powers[degree-2*TileLetters]));
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

    void write_tile(IDEG degree, IDIMN index, IDIMN /*reverse_index*/) noexcept
    {
        assert(0 <= degree && degree <= static_cast<IDEG>(Depth));
        assert(index <= static_cast<IDIMN>(tsi::powers[degree-2*TileLetters]));
        const auto start_of_degree = basis_type::start_of_degree(degree);
        pointer optr = base::out_data + index * tile_width + start_of_degree;
        const_pointer tptr = output_tile.data();
        auto stride = tsi::powers[degree - tile_letters];

        for (DIMN i = 0; i < tile_width; ++i) {
            for (DIMN j = 0; j < tile_width; ++j) {
                optr[i * stride + j] = tptr[i * tile_width + j];
            }
        }

        if (degree < IDEG(Depth)) {
            // Write out reverse data
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
            auto lhs_deg_max = std::min(out_deg-1, old_out_deg);

            if (!assign) {
                update_inplace_lhs<Coeffs>(helper.fwd_write(out_deg),
                                           rhs_unit,
                                           static_cast<IDIMN>(tsi::powers[out_deg]),
                                           op);
            } else if (default_assign && out_deg > rhs_max_deg && lhs_deg_max < lhs_deg_min) {
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
                } else {
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
        OriginalVectors&
        ) const
    {
        if (!lhs.empty() && !rhs.empty()) {
            helper<Coeffs> help(out, lhs, rhs, max_degree);
            fma_impl(help, op, help.out_degree());
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

    template <typename Coeffs, typename Fn>
    LA_INLINE_ALWAYS
    void impl_common(helper_type<Coeffs>& helper,
                     Fn op,
                     IDEG out_deg,
                     IDIMN k,
                     IDIMN k_reverse,
                        IDEG lhs_deg_min,
                        IDEG lhs_deg_max
                        ) const
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
            auto lhs_deg_max = std::min(out_deg-1, helper.lhs_degree());

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
            auto lhs_deg_max = std::min(out_deg-1, old_lhs_deg);

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
        if (max_degree <= 2*tile_info::tile_letters) {
            base::fma(out, lhs, rhs, op, max_degree, orig);
            return;
        }

        if (!lhs.empty() && !rhs.empty()) {
            helper_type<Coeffs> helper(out, lhs, rhs, max_degree);
            fma_impl(helper, op, helper.out_degree());
        }
    }

    template<typename Basis, typename Coeffs, typename Fn, typename OriginalVectors>
    void multiply_inplace(vectors::dense_vector<Basis, Coeffs>& lhs,
                          const vectors::dense_vector<Basis, Coeffs>& rhs,
                          Fn op, DEG max_degree, OriginalVectors& orig) const
    {
        //        std::cout << "BEFORE " << lhs << '\n';
        if (max_degree <= 2*tile_info::tile_letters) {
            base::multiply_inplace(lhs, rhs, op, max_degree, orig);
            return;
        }

        if (!rhs.empty()) {
            helper_type<Coeffs> helper(lhs, rhs, max_degree);
            multiply_inplace_impl(helper, op, helper.out_degree());
        }
        else {
            lhs.clear();
        }
        //        std::cout << "AFTER " << lhs << '\n';
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
class free_tensor : public algebra<
                            free_tensor_basis<n_letters, max_degree>,
                            Coeff,
                            free_tensor_multiplication<n_letters, max_degree>,
                            VectorType,
                            Args...>
{
    typedef free_tensor_multiplication<n_letters, max_degree> multiplication_t;

public:
    /// The basis type.
    typedef free_tensor_basis<n_letters, max_degree> BASIS;
    /// Import of the KEY type.
    typedef typename BASIS::KEY KEY;
    /// The algebra type.
    typedef algebra<BASIS, Coeff, multiplication_t, VectorType, Args...> ALG;

    typedef typename Coeff::SCA SCA;
    typedef typename Coeff::RAT RAT;

    /// The sparse_vector type.
    typedef typename ALG::VECT VECT;

    /// Import of the iterator type.
    typedef typename ALG::iterator iterator;
    /// Import of the constant iterator type.
    typedef typename ALG::const_iterator const_iterator;

public:
    /// Default constructor.
    free_tensor()
    {}

    /// Copy constructor.
    free_tensor(const free_tensor& t)
        : ALG(t)
    {}

    /// Constructs an instance from a shuffle_tensor instance.
    free_tensor(const shuffle_tensor<Coeff, n_letters, max_degree>& t)
    {
        typename shuffle_tensor<Coeff, n_letters, max_degree>::const_iterator i;
        for (i = t.begin(); i != t.end(); ++i) {
            (*this)[i->key()] += i->value();
        }
    }

    /// Constructs an instance from an algebra instance.
    free_tensor(const ALG& a)
        : ALG(a)
    {}

    /// Constructs an instance from a sparse_vector instance.
    free_tensor(const VECT& v)
        : ALG(v)
    {}

    /// Constructs a unidimensional instance from a letter and a scalar.
    free_tensor(LET letter, const SCA& s) : ALG(VECT::basis.keyofletter(letter), s)
    {
    }

    /// Constructor from pointers to data range
    free_tensor(SCA const* begin,
                SCA const* end) : ALG(begin, end)
    {
    }

    /// Explicit unidimensional constructor from a given key (basis element).
    explicit free_tensor(const KEY& k)
        : ALG(k)
    {}

    /// Explicit unidimensional constructor from a given scalar.
    explicit free_tensor(const SCA& s)
        : ALG(VECT::basis.empty_key, s)
    {}

    /// Constructor from pointers to data range with offset
    free_tensor(DIMN
                        offset,
                SCA const* begin, SCA const* end) : ALG(offset, begin, end)
    {
    }

    template <typename B, typename M>
    explicit free_tensor(const algebra<B, Coeff, M, VectorType>& arg)
            : ALG(arg)
    {}

    free_tensor& operator=(const free_tensor&) = default;
    //free_tensor& operator=(free_tensor&&) noexcept = default;

public:
    /// Computes the truncated exponential of a free_tensor instance.
    inline friend free_tensor exp(const free_tensor& arg)
    {
        // Computes the truncated exponential of arg
        // 1 + arg + arg^2/2! + ... + arg^n/n! where n = max_degree
        KEY kunit;
        free_tensor result(kunit);
        for (DEG i = max_degree; i >= 1; --i) {
            result.mul_scal_div(arg, (RAT)i);
            result += (free_tensor)
                    kunit;
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
        KEY kunit;
        typename free_tensor::iterator unit_elt;

        if ((unit_elt = x.find(kunit)) != x.end() && unit_elt->value() != VECT::zero) {
            x.erase(unit_elt);
        }

        for (DEG i = max_degree; i >= 1; --i) {
            result.mul_scal_div(x, static_cast<RAT>(i), max_degree - i + 1);
            result += *this;
        }

        return result;
    }

    /// Inplace version of fmexp
    free_tensor& fmexp_inplace(const free_tensor& exp_arg)
    {
        free_tensor self(*this), x(exp_arg);
        KEY kunit;
        typename free_tensor::iterator unit_elt;

        if ((unit_elt = x.find(kunit)) != x.end() && unit_elt->value() != VECT::zero) {
            x.erase(unit_elt);
        }

        for (DEG i = max_degree; i >= 1; --i) {
            this->mul_scal_div(x, static_cast<RAT>(i), max_degree - i + 1);
            *this += self;
        }

        return *this;
    }

    /// Computes the truncated logarithm of a free_tensor instance.
    inline friend free_tensor log(const free_tensor& arg)
    {
        // Computes the truncated log of arg up to degree max_degree
        // The coef. of the constant term (empty word in the monoid) of arg
        // is forced to 1.
        // log(arg) = log(1+x) = x - x^2/2 + ... + (-1)^(n+1) x^n/n.
        // max_degree must be > 0
        KEY kunit;
        free_tensor tunit(kunit);
        free_tensor x(arg);
        iterator it = x.find(kunit);
        if (it != x.end()) {
            x.erase(it);
        }
        free_tensor result;

        for (DEG i = max_degree; i >= 1; --i) {
            if (i % 2 == 0) {
                result.sub_scal_div(tunit, (RAT)i);
            }
            else {
                result.add_scal_div(tunit, (RAT)i);
            }
            result *= x;
        }

        return result;
    }

    /// Computes the truncated inverse of a free_tensor instance.
    inline friend free_tensor inverse(const free_tensor& arg)
    {
        // Computes the truncated inverse of arg up to degree max_degree
        // An exception is thrown if the leading term is zero.
        // the module assumes
        // (a+x)^(-1) = (a(1+x/a))^(-1)
        //  = a^(-1)(1 - x/a + x^2/a^2 + ... + (-1)^(n) x^n/a^n)
        // = a^(-1) - x/a*[a^(-1)(1 - x/a + x^2/a^2 + ... + (-1)^(n)
        // x^(n-1)/a^(n-1)))]. S_n = a^(-1) + z S_{n-1}; z = - x/a ; S_0 = a^(-1)
        // max_degree must be > 0

        static KEY kunit;
        SCA a(0);
        free_tensor x, z(a);

        const_iterator it(arg.find(kunit));
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
        free_tensor free_tensor_a_inverse(SCA(1) / a), result(free_tensor_a_inverse);
        // z := - x/a
        z.sub_scal_div(x, a);
        // the iteration
        for (DEG i = 0; i != max_degree; ++i) {
            result = free_tensor_a_inverse + z * result;
        }
        return result;
    }

    /// Computes the truncated inverse of a free_tensor instance.
    inline friend free_tensor reflect(const free_tensor& arg)
    {
        // Computes the alternating reflection of arg up to degree max_degree
        // For group-like elements this is the same as the inverse
        free_tensor ans(SCA(0));
        for (const_iterator it = arg.begin(); it != arg.end(); ++it) {
            KEY old_key = it->key();
            SCA old_value = it->value();
            ans[old_key.reverse()] = (old_key.size() % 2) ? SCA(0) - old_value : old_value;
        }
        return ans;
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
        ar& boost::serialization::base_object<ALG>(*this);
    }
#endif
};

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
    shuffle_tensor(const free_tensor<Coeff, n_letters, max_degree>& t)
    {
        typename free_tensor<Coeff, n_letters, max_degree>::const_iterator i;
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

}// namespace alg
// Include once wrapper
#endif// DJC_COROPA_LIBALGEBRA_TENSORH_SEEN

// EOF.
