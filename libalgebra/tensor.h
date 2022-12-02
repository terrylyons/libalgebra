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

//#include <omp.h>

#include <algorithm>
#include <unordered_set>

#include <boost/align/aligned_alloc.hpp>
#include <boost/align/aligned_allocator.hpp>
#include <boost/call_traits.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/functional/hash.hpp>

#include "algebra.h"
#include "area_tensor_basis.h"
#include "base_vector.h"
#include "dense_storage.h"
#include "dense_vector.h"
#include "detail/integer_maths.h"
#include "detail/level_walkers.h"
#include "detail/reversing_permutation.h"
#include "half_shuffle_tensor_basis.h"
#include "tensor_basis.h"

#define LA_RESTRICT __restrict
#define LA_INLINE_ALWAYS __attribute__((always_inline))

#define LA_ALIGNAS(BYTES) alignas(BYTES)
#ifndef LA_CACHELINE_BYTES
#define LA_CACHELINE_BYTES 64
#endif

#ifndef LIBALGEBRA_L1_CACHE_SIZE
#define LIBALGEBRA_L1_CACHE_SIZE 32768// 32KB should be fairly standard
#endif

#define LA_DEFAULT_TILE_PARAM(WIDTH, SCALAR) ::alg::dtl::tensor_tile_letters_helper<WIDTH, LIBALGEBRA_L1_CACHE_SIZE / (2 * sizeof(SCALAR))>::num_letters

namespace alg {

namespace dtl {

struct default_signer {
    explicit default_signer(DEG degree) : sign(degree % 2 == 0 ? 1 : -1)
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

    explicit non_signing_signer(DEG deg) {}
};

template<DEG Width, DIMN TargetSize, bool WidthSquareFits = (Width * Width < TargetSize)>
struct tensor_tile_letters_helper {
    static constexpr IDEG log_target = IDEG(::alg::integer_maths::logN(TargetSize, Width));
#ifdef LIBALGEBRA_MAX_TILE_LETTERS
    static constexpr IDEG num_letters = std::min(
            IDEG(LIBALGEBRA_MAX_TILE_LETTERS),
            (log_target >= 2) ? log_target / 2 : -1);
#else
    static constexpr IDEG num_letters = (log_target >= 2) ? log_target / 2 : -1;
#endif
};

template<DEG Width, DIMN TargetSize>
struct tensor_tile_letters_helper<Width, TargetSize, false> {
    /*
     * If we're here, it is because Width^2 >= TargetSize.
     * In this case we're looking for n such that (Width/2^n)^2 < TargetSize,
     * so n = -log_2(TargetSize/Width^2)/2, but this won't work since TargetSize/Width^2
     * (in integer maths) will be 0.
     * Instead, find m with m >= log_2(Width^2/TargetSize)/2, then take n=-m
     */
    static constexpr IDEG log_target = ::alg::integer_maths::logN(IDIMN(Width * Width) / IDIMN(TargetSize), IDEG(2));
    static constexpr IDEG num_letters = -((log_target % 2 == 0) ? log_target / 2 : 1 + log_target / 2);
};

template<typename S, DIMN Size>
struct data_tile {
    static_assert(Size > 0, "size must be non-zero");
    static constexpr DIMN size = Size;
    S* data;

    data_tile()
    {
        data = static_cast<S*>(boost::alignment::aligned_alloc(LA_CACHELINE_BYTES, Size * sizeof(S)));
        if (!std::is_trivially_default_constructible<S>::value) {
            auto* ptr = data;
            for (DIMN i = 0; i < Size; ++i, ++ptr) {
                ::new (ptr) S();
            }
        }
    }

    ~data_tile()
    {
        boost::alignment::aligned_free(data);
    }
};

template<DEG Width, DEG Depth, IDEG TileLetters = 0>
struct tile_details {
    static constexpr IDEG tile_letters = (TileLetters > 0) ? TileLetters : 1;
    static constexpr IDIMN tile_width = (TileLetters >= 0)
            ? integer_maths::power(IDIMN(Width), tile_letters)
            : IDIMN(Width) / integer_maths::power(IDIMN(2), -TileLetters);
    static constexpr IDIMN num_subtiles = (TileLetters >= 0)
            ? 1
            : (integer_maths::power(IDIMN(2), -TileLetters)
               + (Width % integer_maths::power(IDIMN(2), -TileLetters) == 0 ? 0 : 1));

    static constexpr IDIMN tile_stride = (TileLetters >= 0) ? tile_width : Width;
    static constexpr IDIMN tile_size = tile_width * tile_width;
    static constexpr IDIMN tile_shift = integer_maths::power(Width, tile_letters - 1);
};

template<DEG Width, DEG Depth, typename Coeffs, IDEG TileLetters>
class central_tile_helper : public tile_details<Width, Depth, TileLetters>
{
public:
    using scalar_type = typename Coeffs::S;
    using pointer = scalar_type*;
    using const_pointer = const scalar_type*;
    using size_type = std::size_t;
    using index_type = std::ptrdiff_t;
    using degree_type = IDEG;
    using tile_info = tile_details<Width, Depth, TileLetters>;
    using tsi = tensor_size_info<Width>;

private:
    static constexpr index_type tile_width = static_cast<index_type>(tile_info::tile_width);
    data_tile<scalar_type, tile_info::tile_size> tile;

public:
    constexpr static IDIMN pointer_offset(IDEG degree, IDIMN index_word, IDIMN subtile_i, IDIMN subtile_j) noexcept
    {
        return (degree == 0) ? 0 : tsi::degree_sizes[degree - 1] + index_word * tile_info::tile_stride + (subtile_i * tsi::powers[degree - tile_info::tile_letters] + subtile_j) * tile_info::tile_width;
    }

    static bool boundary_subtile(IDIMN index) noexcept
    {
        return (tile_info::num_subtiles != 1) && (Width % tile_width != 0) && (DIMN(index) == tile_info::num_subtiles - 1);
    }

    static index_type subtile_bound(index_type subtile) noexcept
    {
        return (boundary_subtile(subtile) ? Width % tile_width : tile_width);
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
    }

protected:
    void read_tile_impl(const_pointer in_p, index_type lhs_stride, index_type ibound, index_type jbound) noexcept
    {
        for (index_type i = 0; i < ibound; ++i) {
            for (index_type j = 0; j < jbound; ++j) {
                tile.data[i * tile_width + j] = in_p[i * lhs_stride + j];
            }
        }
    }

    void write_tile_impl(pointer out_p, index_type lhs_stride, index_type ibound, index_type jbound) const noexcept
    {
        for (index_type i = 0; i < ibound; ++i) {
            for (index_type j = 0; j < jbound; ++j) {
                out_p[i * lhs_stride + j] = tile.data[i * tile_width + j];
            }
        }
    }

    void write_tile_reverse_impl(pointer out_p, index_type lhs_stride, index_type ibound, index_type jbound) const noexcept
    {
        using perm = reversing_permutation<Width, tile_info::tile_letters>;
        for (index_type i = 0; i < ibound; ++i) {
            for (index_type j = 0; j < jbound; ++j) {
                out_p[i * lhs_stride + j] = tile.data[perm::permute_idx(j) * tile_width + perm::permute_idx(i)];
            }
        }
    }

private:
    static std::vector<int, boost::alignment::aligned_allocator<int, LA_CACHELINE_BYTES>> setup_reverser()
    {
        using perm = reversing_permutation<Width, tile_info::tile_letters>;
        std::vector<int, boost::alignment::aligned_allocator<int, LA_CACHELINE_BYTES>> result;
        result.reserve(tile_width);
        for (DIMN i = 0; i < tile_width; ++i) {
            result.push_back((int)perm::permute_idx(i));
        }
        return result;
    }

public:
    static const int* reverser()
    {
        static const auto reverse = setup_reverser();
        return reverse.data();
    }

    pointer tile_ptr() noexcept
    {
        return static_cast<pointer>(tile.data);
    }
    const_pointer tile_ptr() const noexcept
    {
        return static_cast<const_pointer>(tile.data);
    }

    void read_tile(const_pointer src_p,
                   degree_type degree,
                   index_type index,
                   index_type subtile_i = 0,
                   index_type subtile_j = 0) noexcept
    {
        const auto stride = tsi::powers[degree - tile_info::tile_letters];
        const_pointer in_p = src_p + pointer_offset(degree, index, subtile_i, subtile_j);
        read_tile_impl(in_p, stride, subtile_bound(subtile_i), subtile_bound(subtile_j));
    }

    void write_tile(pointer dst_p,
                    degree_type degree,
                    index_type index,
                    index_type subtile_i = 0,
                    index_type subtile_j = 0) const noexcept
    {
        const auto stride = tsi::powers[degree - tile_info::tile_letters];
        pointer out_p = dst_p + pointer_offset(degree, index, subtile_i, subtile_j);
        write_tile_impl(out_p, stride, subtile_bound(subtile_i), subtile_bound(subtile_j));
    }

    void write_tile_reverse(pointer dst_p,
                            degree_type degree,
                            index_type index,
                            index_type subtile_i = 0,
                            index_type subtile_j = 0) const noexcept
    {
        const auto stride = tsi::powers[degree - tile_info::tile_letters];
        pointer out_p = dst_p + pointer_offset(degree, index, subtile_j, subtile_i);
        write_tile_reverse_impl(out_p, stride, subtile_bound(subtile_i), subtile_bound(subtile_j));
    }

    void reset_tile() noexcept
    {
        std::fill(tile.data, tile.data + tile_info::tile_size, Coeffs::zero);
    }
};

template<DEG Width, DEG Depth, typename Coeffs, IDEG TileLetters>
class tiled_antipode_helper : public central_tile_helper<Width, Depth, Coeffs, TileLetters>
{
    using base = central_tile_helper<Width, Depth, Coeffs, TileLetters>;
    using tile_info = tile_details<Width, Depth, TileLetters>;

public:
    using typename base::const_pointer;
    using typename base::pointer;

private:
    const_pointer src_ptr;
    const_pointer src_reverse_ptr = nullptr;
    pointer dst_ptr;
    pointer dst_reverse_ptr = nullptr;

    template<typename B>
    using dense_tensor = ::alg::vectors::dense_vector<B, Coeffs>;

    using dense_free_tensor = ::alg::vectors::dense_vector<free_tensor_basis<Width, Depth>, Coeffs>;

public:
    using base::boundary_subtile;
    using base::num_subtiles;
    using base::pointer_offset;
    using base::subtile_bound;
    using base::tile_ptr;
    using tile_info::tile_letters;
    using tile_info::tile_size;
    using tile_info::tile_width;

    void read_tile(IDEG degree, IDIMN index, IDIMN subtile_i, IDIMN subtile_j)
    {
        base::read_tile(src_ptr, degree, index, subtile_i, subtile_j);
    }
    void write_tile(IDEG degree, IDIMN index, IDIMN subtile_i, IDIMN subtile_j) const
    {
        base::write_tile(dst_ptr, degree, index, subtile_i, subtile_j);
    }

    const_pointer reverse_src_ptr() const noexcept { return src_reverse_ptr; }
    pointer reverse_dst_ptr() const noexcept { return dst_reverse_ptr; }

    tiled_antipode_helper(pointer dst_p, const_pointer src_p)
        : src_ptr(src_p), dst_ptr(dst_p)
    {}

    template<typename B>
    tiled_antipode_helper(dense_tensor<B>& result, const dense_tensor<B>& arg)
    {
        src_ptr = arg.as_ptr();
        dst_ptr = result.as_mut_ptr();
    }

    tiled_antipode_helper(dense_free_tensor& result, const dense_free_tensor& arg)
    {
        src_ptr = arg.as_ptr();
        dst_ptr = result.as_mut_ptr();
    }
};

template<DEG Width, DEG MaxDepth, typename Coeffs, typename Signer,
         IDEG TileLetters = LA_DEFAULT_TILE_PARAM(Width, typename Coeffs::S)>
class tiled_inverse_operator
{
    using scalar_type = typename Coeffs::S;

    using BASIS = free_tensor_basis<Width, MaxDepth>;
    using tsi = dtl::tensor_size_info<Width>;
    using tile_info = tile_details<Width, MaxDepth, TileLetters>;

    template<typename C>
    using helper_type = tiled_antipode_helper<Width, MaxDepth, C, TileLetters>;

    template<DEG Degree>
    struct untiled_compute {
        template<typename S>
        static void eval(S* LA_RESTRICT dst_ptr, const S* LA_RESTRICT src_ptr, DEG current_degree)
        {
            if (current_degree >= Degree) {
                using perm = dtl::reversing_permutation<Width, Degree>;

                constexpr DIMN bound = ::alg::integer_maths::power(DIMN(Width), Degree);
                Signer signer(Degree);

                constexpr DIMN offset = (Degree == 0) ? 0 : integer_maths::sum_powers(DIMN(Width), Degree - 1);
                scalar_type* LA_RESTRICT out_ptr = dst_ptr + offset;
                const scalar_type* LA_RESTRICT in_ptr = src_ptr + offset;
                for (DIMN i = 0; i < bound; ++i) {
                    out_ptr[perm::permute_idx(i)] = signer(in_ptr[i]);
                }
            }
        }
    };

    template<typename C>
    using pointer = typename C::S*;
    template<typename C>
    using const_pointer = const typename C::S*;

public:
    static constexpr DEG block_letters = tile_info::tile_letters;

    static void untiled_cases(scalar_type* LA_RESTRICT dst_ptr, const scalar_type* LA_RESTRICT src_ptr, DEG current_degree) noexcept
    {
        dtl::increasing_level_walker<untiled_compute, 2 * tile_info::tile_letters>::eval(dst_ptr, src_ptr, current_degree);
    }

    static void sign_and_permute(scalar_type* LA_RESTRICT tile, Signer& signer) noexcept
    {
        for (DIMN i = 0; i < tile_info::tile_size; ++i) {
            tile[i] = signer(tile[i]);
        }

        using perm = reversing_permutation<Width, 2 * block_letters>;
        for (DIMN i = 0; i < tile_info::tile_size; ++i) {
            auto j = perm::permute_idx(i);
            if (j > i) {
                std::swap(tile[j], tile[i]);
            }
        }
    }

    static void permute_level_tiled(helper_type<Coeffs>& helper, DEG out_deg) noexcept
    {
        const auto mid_deg = out_deg - 2 * tile_info::tile_letters;
        Signer signer(out_deg);

        for (IDIMN middle_index = 0; middle_index < IDIMN(tsi::powers[mid_deg]); ++middle_index) {
            const auto reverse_middle_index = helper.reverse_key(mid_deg, middle_index);

            for (IDIMN subtile_i = 0; subtile_i < tile_info::num_subtiles; ++subtile_i) {
                for (IDIMN subtile_j = 0; subtile_j < tile_info::num_subtiles; ++subtile_j) {
                    helper.read_tile(out_deg, middle_index, subtile_i, subtile_j);
                    sign_and_permute(helper.tile_ptr(), signer);
                    helper.write_tile(out_deg, reverse_middle_index, subtile_i, subtile_j);
                }
            }
        }
    }

    static void tile_cases(helper_type<Coeffs>& helper, DEG current_degree) noexcept
    {
        for (DEG out_deg = 2 * tile_info::tile_letters + 1; out_deg <= current_degree; ++out_deg) {
            permute_level_tiled(helper, out_deg);
        }
    }

    void operator()(const scalar_type* src_ptr, scalar_type* dst_ptr, const DEG curr_degree) const noexcept
    {

        if (src_ptr == nullptr)// if pointer to source is null
        {
            return;
        }
        if (dst_ptr == nullptr) {
            return;
        }

        untiled_cases(dst_ptr, src_ptr, curr_degree);
        if (curr_degree > 2 * tile_info::tile_letters) {
            helper_type<Coeffs> helper(dst_ptr, src_ptr);
            tile_cases(helper, curr_degree);
        }
    }

    template<typename C = Coeffs>
    static void apply_pointers(const_pointer<C> src_ptr, pointer<C> dst_ptr, DEG max_degree)
    {
        if (src_ptr == nullptr)// if pointer to source is null
        {
            return;
        }
        if (dst_ptr == nullptr) {
            return;
        }

        untiled_cases(dst_ptr, src_ptr, max_degree);
        if (max_degree > 2 * tile_info::tile_letters) {
            helper_type<Coeffs> helper(dst_ptr, src_ptr);
            tile_cases(helper, max_degree);
        }
    }

    template<typename Vector>
    static void apply(const Vector& src, Vector& result, DEG max_degree = MaxDepth)
    {
        using ring = typename Vector::coefficient_ring;
        for (auto item : src) {
            auto key = item.key();
            auto deg = key.size();
            if (max_degree == MaxDepth && deg <= max_degree) {
                if (key.size() % 2 == 0) {
                    result[key.reverse()] = item.value();
                }
                else {
                    result[key.reverse()] = ring::uminus(item.value());
                }
            }
        }
    }

    template<typename B, typename C>
    static void apply(const vectors::dense_vector<B, C>& src, vectors::dense_vector<B, C>& dst, DEG max_depth = MaxDepth)
    {
        if (src.dimension() == 0) {
            return;
        }
        assert(max_depth <= MaxDepth);
        dst.resize_to_degree(std::min(src.degree(), max_depth));
        untiled_cases(dst.as_mut_ptr(), src.as_ptr(), std::min(src.degree(), max_depth));
        if (max_depth > 2 * helper_type<C>::tile_letters) {
            helper_type<C> helper(dst, src);
            tile_cases(helper, std::min(src.degree(), max_depth));
        }
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

    using index_type = typename central_tile_helper<Width, Depth, Coeffs, 0>::index_type;

    template<typename B>
    using dense_tensor = vectors::dense_vector<B, coefficient_ring>;

    std::vector<const_pointer> m_lhs_levels;
    std::vector<const_pointer> m_rhs_levels;
    std::vector<pointer> m_out_levels;
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

protected:
    template<typename P>
    static void setup_level_list(std::vector<P>& levels, IDEG depth)
    {
        for (IDEG i = 1; i <= depth; ++i) {
            levels.push_back(levels.back() + tsi::powers[i - 1]);
        }
    }

private:
    void setup_levels()
    {
        m_lhs_levels.reserve(lhs_deg + 1);
        m_rhs_levels.reserve(rhs_deg + 1);
        m_out_levels.reserve(out_deg + 1);

        m_lhs_levels.push_back(lhs_data == nullptr ? out_data : lhs_data);
        setup_level_list(m_lhs_levels, lhs_deg);
        m_rhs_levels.push_back(rhs_data);
        setup_level_list(m_rhs_levels, rhs_deg);
        m_out_levels.push_back(out_data);
        setup_level_list(m_out_levels, out_deg);
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
        setup_levels();
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
        setup_levels();
    }

    constexpr IDEG lhs_degree() const noexcept { return lhs_deg; }
    constexpr IDEG rhs_degree() const noexcept { return rhs_deg; }
    constexpr IDEG out_degree() const noexcept { return out_deg; }
    constexpr bool is_inplace() const noexcept { return lhs_data == nullptr; }

    reference out_unit() noexcept { return *m_out_levels[0]; }
    const_reference left_unit() const noexcept { return *m_lhs_levels[0]; }
    const_reference right_unit() const noexcept { return *m_rhs_levels[0]; }
    const_pointer left_fwd_read(IDEG d, IDIMN offset = 0) const noexcept
    {
        assert(d >= 0 && d <= lhs_deg);
        assert(offset + basis_type::start_of_degree(d) <= basis_type::start_of_degree(d + 1));
        //        return ((lhs_data == nullptr) ? out_data : lhs_data) + offset + basis_type::start_of_degree(d);
        assert(offset < IDIMN(tsi::powers[d]));
        return m_lhs_levels[d] + offset;
    }
    const_pointer right_fwd_read(IDEG d, IDIMN offset = 0) const noexcept
    {
        assert(d >= 0 && d <= rhs_deg);
        assert(offset + basis_type::start_of_degree(d) <= basis_type::start_of_degree(d + 1));
        //        return rhs_data + offset + basis_type::start_of_degree(d);
        assert(offset < IDIMN(tsi::powers[d]));
        return m_rhs_levels[d] + offset;
    }
    pointer fwd_write(IDEG d) const noexcept
    {
        assert(d >= 0 && d <= out_deg);
        //        return out_data + basis_type::start_of_degree(d);
        return m_out_levels[d];
    }

    std::pair<DIMN, DIMN> range_size(IDEG lhs, IDEG rhs) const noexcept
    {
        return std::pair<DIMN, DIMN>(tsi::powers[lhs], tsi::powers[rhs]);
    }
};

template<DEG Width, DEG Depth, typename Coeffs, IDEG TileLetters>
class tiled_free_tensor_multiplication_helper
    : public free_tensor_multiplication_helper<Width, Depth, Coeffs>,
      public central_tile_helper<Width, Depth, Coeffs, TileLetters>
{
    using base = free_tensor_multiplication_helper<Width, Depth, Coeffs>;
    using tile_info = tile_details<Width, Depth, TileLetters>;
    using base_helper = central_tile_helper<Width, Depth, Coeffs, TileLetters>;
    using tsi = tensor_size_info<Width>;

public:
    using scalar_type = typename Coeffs::SCA;
    using pointer = scalar_type*;
    using const_pointer = const scalar_type*;
    using reference = scalar_type&;
    using const_reference = const scalar_type&;
    using typename base_helper::index_type;

private:
    template<typename B>
    using dense_tensor = vectors::dense_vector<B, Coeffs>;

    using basis_type = tensor_basis<Width, Depth>;

    std::vector<const_pointer> m_lhs_rev_levels;
    std::vector<pointer> m_out_rev_levels;
    //    std::vector<scalar_type> left_read_tile;
    dtl::data_tile<scalar_type, tile_info::tile_width> left_read_tile;
    //    std::vector<scalar_type> right_read_tile;
    dtl::data_tile<scalar_type, tile_info::tile_width> right_read_tile;
    //    std::vector<scalar_type> output_tile;
    std::vector<scalar_type> reverse_data;
    const_pointer left_reverse_read_ptr = nullptr;
    pointer reverse_write_ptr = nullptr;

    template<typename B>
    static void fill_reverse_data(pointer out, const dense_tensor<B>& lhs, IDEG degree)
    {
        assert(0 < degree && DEG(degree) <= Depth);
        dtl::tiled_inverse_operator<Width, (Depth > 0) ? Depth - 1 : 0, Coeffs, dtl::non_signing_signer, TileLetters> reverser;
        reverser(lhs.as_ptr(), out, DEG(degree - 1));
    }

    template<typename B>
    void allocate_and_fill_reverse_data(const dense_tensor<B>& lhs, IDEG degree)
    {
        reverse_data.resize(basis_type::start_of_degree(DEG(degree)));
        fill_reverse_data(reverse_data.data(), lhs, degree);
        left_reverse_read_ptr = reverse_data.data();
    }

    template<typename B>
    void setup_reverse_read(const dense_tensor<B>& lhs)
    {
        if (base::lhs_deg > 0) {
            allocate_and_fill_reverse_data(lhs, base::lhs_deg);
        }
    }

    void setup_reverse_read(
            const dense_tensor<free_tensor_basis<Width, Depth>>& lhs)
    {
        if (base::lhs_deg > 0) {
            const auto& vect_reverse_data = lhs.reverse_data();
            if (vect_reverse_data) {
                assert(vect_reverse_data.size() == basis_type::start_of_degree(DEG(base::lhs_deg)));
                left_reverse_read_ptr = vect_reverse_data.begin();
            }
            else {
                allocate_and_fill_reverse_data(lhs, base::lhs_deg);
            }
        }
    }

    template<typename B>
    void setup_reverse_write(dense_tensor<B>& out)
    {}

    void setup_reverse_write_impl(dense_tensor<free_tensor_basis<Width, Depth>>& out, DEG degree)
    {
        auto& vect_reverse_data = out.reverse_data();
        auto size = basis_type::start_of_degree(DEG(degree));
        if (size > vect_reverse_data.size()) {
            vect_reverse_data.resize(size);
        }
        reverse_write_ptr = vect_reverse_data.begin();
        vect_reverse_data.validate();
    }

    void setup_reverse_write(dense_tensor<free_tensor_basis<Width, Depth>>& out)
    {
        if (base::out_deg > 0) {
            setup_reverse_write_impl(out, base::out_deg);
        }
    }

    template<typename B>
    void setup_reverse_readwrite(dense_tensor<B>& out)
    {
        if (base::lhs_deg > 0) {
            allocate_and_fill_reverse_data(out, base::lhs_deg);
        }
    }

    void setup_reverse_readwrite(dense_tensor<free_tensor_basis<Width, Depth>>& out)
    {
        // if the old lhs_deg is larger than then we don't necessarily need to do
        // any resizing, but we do need to make sure the pointers are set up
        // correctly, and that the reverse data is valid.
        const auto deg = std::max(base::out_deg, base::lhs_deg);
        if (deg > 0) {
            auto& vect_reverse_data = out.reverse_data();
            bool was_valid = static_cast<bool>(vect_reverse_data);
            setup_reverse_write_impl(out, deg);

            // Rather than allocate a separate buffer for reverse data, we
            // can just use the reverse data buffer on the out vector.
            if (!was_valid && base::lhs_deg > 0) {
                fill_reverse_data(reverse_write_ptr, out, base::lhs_deg);
            }
            left_reverse_read_ptr = reverse_write_ptr;
        }
    }

    void setup_reverse_levels()
    {
        m_lhs_rev_levels.reserve(base::lhs_deg);
        m_out_rev_levels.reserve(base::out_deg);

        m_lhs_rev_levels.push_back(left_reverse_read_ptr);
        base::setup_level_list(m_lhs_rev_levels, (base::lhs_deg > 0) ? base::lhs_deg - 1 : 0);
        if (reverse_write_ptr != nullptr) {
            m_out_rev_levels.push_back(reverse_write_ptr);
            base::setup_level_list(m_out_rev_levels, (base::out_deg > 0) ? base::out_deg - 1 : 0);
        }
    }

public:
    using base_helper::boundary_subtile;
    using base_helper::pointer_offset;
    using base_helper::reverse_key;
    using base_helper::reverser;

    using tile_info::tile_letters;
    using tile_info::tile_shift;
    using tile_info::tile_size;
    using tile_info::tile_width;

    template<typename B1, typename B2, typename B3>
    tiled_free_tensor_multiplication_helper(dense_tensor<B1>& out,
                                            const dense_tensor<B2>& lhs,
                                            const dense_tensor<B3>& rhs,
                                            DEG max_degree)
        : base(out, lhs, rhs, max_degree)
    {
        setup_reverse_read(lhs);
        setup_reverse_write(out);
        setup_reverse_levels();
    }

    template<typename B1, typename B2>
    tiled_free_tensor_multiplication_helper(dense_tensor<B1>& lhs, const dense_tensor<B2>& rhs, DEG max_degree)
        : base(lhs, rhs, max_degree)
    {
        setup_reverse_readwrite(lhs);
        setup_reverse_levels();
    }

    pointer out_tile_ptr() noexcept
    {
        return base_helper::tile_ptr();
    }
    const_pointer left_read_tile_ptr() const noexcept
    {
        return left_read_tile.data;
    }
    const_pointer right_read_tile_ptr() const noexcept
    {
        return right_read_tile.data;
    }

    void read_left_tile(IDEG degree, IDIMN index, IDIMN subtile_i = 0) noexcept
    {
        const auto* ptr_begin = left_reverse_read_ptr + pointer_offset(degree, index, 0, subtile_i);
        if (boundary_subtile(subtile_i)) {
            const auto mid = Width % tile_width;
            std::copy(ptr_begin, ptr_begin + mid, left_read_tile.data);
            std::fill(left_read_tile.data + mid, left_read_tile.data + tile_width, Coeffs::zero);
        }
        else {
            std::copy(ptr_begin, ptr_begin + tile_width, left_read_tile.data);
        }
    }

    void read_left_tile(const_pointer src, IDIMN count = tile_width) noexcept
    {
        if (count == tile_width) {
            std::copy(src, src + count, left_read_tile.data);
        }
        else {
            std::copy(src, src + count, left_read_tile.data);
            std::fill(left_read_tile.data + count, left_read_tile.data + tile_width, Coeffs::zero);
        }
    }
    void read_right_tile(IDEG degree, IDIMN index, IDIMN subtile_j = 0) noexcept
    {
        const auto* ptr_begin = base::rhs_data + pointer_offset(degree, index, 0, subtile_j);
        if (boundary_subtile(subtile_j)) {
            const auto mid = Width % tile_width;
            std::copy(ptr_begin, ptr_begin + mid, right_read_tile.data);
            std::fill(right_read_tile.data + mid, right_read_tile.data + tile_width, Coeffs::zero);
        }
        else {
            std::copy(ptr_begin, ptr_begin + tile_width, right_read_tile.data);
        }
    }
    void read_right_tile(const_pointer src, IDIMN count = tile_width) noexcept
    {
        // TODO: make resilient to subtiles
        if (count == tile_width) {
            std::copy(src, src + count, right_read_tile.data);
        }
        else {
            std::copy(src, src + count, right_read_tile.data);
            std::fill(right_read_tile.data + count, right_read_tile.data + tile_width, Coeffs::zero);
        }
    }

    void reset_tile(IDEG degree, IDIMN index, IDIMN /*reverse_index*/, IDIMN subtile_i = 0, IDIMN subtile_j = 0) noexcept
    {
#if 0
        if (base::lhs_data != nullptr) {
            assert(0 <= degree && degree <= static_cast<IDEG>(Depth));
            assert(index < static_cast<IDEG>(tsi::powers[degree - 2 * tile_letters]));
            base_helper::read_tile(base::out_data, degree, index, subtile_i, subtile_j);
        }
        else {
#endif
        base_helper::reset_tile();
#if 0
        }
#endif
    }

    void write_tile(IDEG degree, IDIMN index, IDIMN reverse_index, IDIMN subtile_i = 0, IDIMN subtile_j = 0) noexcept
    {
        assert(0 <= degree && degree <= static_cast<IDEG>(Depth));
        assert(index <= static_cast<IDIMN>(tsi::powers[degree - 2 * tile_letters]));
        assert(reverse_index <= static_cast<IDIMN>(tsi::powers[degree - 2 * tile_letters]));
        base_helper::write_tile(base::out_data, degree, index, subtile_i, subtile_j);

        if (reverse_write_ptr != nullptr && degree < base::out_deg) {
            // Write out reverse data
            base_helper::write_tile_reverse(reverse_write_ptr, degree, reverse_index, subtile_i, subtile_j);
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

    const_pointer left_fwd_read_ptr(IDEG degree, IDIMN index, IDIMN subtile_i, IDIMN subtile_j = 0) const noexcept
    {
        return base::left_fwd_read(degree, index * tile_info::tile_stride + subtile_i * tile_info::tile_width * tsi::powers[degree - tile_letters]);
    }
    const_pointer left_rev_read_ptr(IDEG degree, IDIMN index, IDIMN subtile_i, IDIMN subtile_j = 0) const noexcept
    {
        return m_lhs_rev_levels[degree] + index * tile_info::tile_stride + subtile_i * tile_width * tsi::powers[degree - tile_letters];
    }
    const_pointer right_fwd_read_ptr(IDEG degree, IDIMN index, IDIMN subtile_j, IDIMN subtile_i = 0) const noexcept
    {
        return base::right_fwd_read(degree, index * tile_info::tile_stride + subtile_j * tile_info::tile_width);
    }
};

}// namespace dtl

template<DEG Width, DEG Depth>
class free_tensor_multiplier
{
public:
    using basis_type = tensor_basis<Width, Depth>;
    using key_type = typename basis_type::KEY;
    using index_type = std::ptrdiff_t;
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
    using index_type = std::ptrdiff_t;

    using basis_type = tensor_basis<Width, Depth>;

protected:
    template<typename B, typename Coeffs>
    static void update_reverse_data(vectors::dense_vector<B, Coeffs>& out, DEG max_degree)
    {}

    template<typename Coeffs>
    static void update_reverse_data(vectors::dense_vector<free_tensor_basis<Width, Depth>, Coeffs>& out, DEG max_degree)
    {
        out.construct_reverse_data(max_degree);
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
        if (lhs.dimension() == 0 && rhs.dimension() != 0) {
            helper<Coeffs> help(out, lhs, rhs, max_degree);
            fma_impl(help, op, help.out_degree());
            update_reverse_data(out, help.out_degree());
        }
    }

    template<typename Basis, typename Coeffs, typename Fn, typename OriginalVectors>
    std::enable_if_t<std::is_base_of<basis_type, Basis>::value>
    multiply_inplace(vectors::dense_vector<Basis, Coeffs>& lhs,
                     const vectors::dense_vector<Basis, Coeffs>& rhs,
                     Fn op,
                     DEG max_degree, OriginalVectors&) const
    {
        if (rhs.dimension() != 0) {
            helper<Coeffs> help(lhs, rhs, max_degree);
            multiply_inplace_impl(help, op, help.out_degree());
            update_reverse_data(lhs, help.out_degree());
        }
        else {
            lhs.clear();
        }
    }
};

template<DEG Width, DEG Depth, IDEG TileLetters = 0>
class tiled_free_tensor_multiplication
    : public traditional_free_tensor_multiplication<Width, Depth>
{
    using base = traditional_free_tensor_multiplication<Width, Depth>;

    bool debug = false;

    template<typename C>
    using helper_type = dtl::tiled_free_tensor_multiplication_helper<
            Width, Depth, C, TileLetters == 0 ? LA_DEFAULT_TILE_PARAM(Width, typename C::S) : TileLetters>;

    using tsi = dtl::tensor_size_info<Width>;

    template<typename C>
    using pointer = typename C::S* LA_RESTRICT;

    template<typename C>
    using const_pointer = const typename C::S* LA_RESTRICT;

    template<typename C>
    using const_reference = const typename C::S&;

public:
    using index_type = std::ptrdiff_t;

private:
    template<typename C, typename Fn>
    LA_INLINE_ALWAYS static void linear_fwd_zipped_mul(pointer<C> out_p,
                                                       const_pointer<C> lsrc,
                                                       const_pointer<C> rsrc,
                                                       const_reference<C> lunit,
                                                       const_reference<C> runit,
                                                       index_type bound,
                                                       Fn op) noexcept
    {
        for (index_type i = 0; i < bound; ++i) {
            out_p[i] += op(lsrc[i] * runit) + op(lunit * rsrc[i]);
        }
    }

    template<typename C, typename Fn>
    LA_INLINE_ALWAYS static void
    linear_fwd_zipped_mul_inplace(
            pointer<C> out_p,
            const_pointer<C> lsrc,
            const_pointer<C> rsrc,
            const_reference<C> lunit,
            const_reference<C> runit,
            index_type bound,
            Fn op) noexcept
    {
        for (index_type i = 0; i < bound; ++i) {
            out_p[i] = op(lsrc[i] * runit) + op(lunit * rsrc[i]);
        }
    }

    template <typename C, typename Fn>
    LA_INLINE_ALWAYS static void
    linear_fwd_mul_inplace_left(pointer<C> optr, const_pointer<C> lptr, const_reference<C> runit, index_type bound, Fn op) noexcept
    {
        for (index_type i=0; i<bound; ++i) {
            optr[i] = op(lptr[i]*runit);
        }
    }

    template <typename C, typename Fn>
    LA_INLINE_ALWAYS static void
    linear_fwd_mul_left(pointer<C> optr, const_pointer<C> lptr, const_reference<C> runit, index_type bound, Fn op) noexcept
    {
        for (index_type i=0; i<bound; ++i) {
            optr[i] += op(lptr[i]*runit);
        }
    }

    template <typename C, typename Fn>
    LA_INLINE_ALWAYS static void
    linear_fwd_mul_inplace_right(pointer<C> optr, const_pointer<C> rptr, const_reference<C> lunit, index_type bound, Fn op) noexcept
    {
        for (index_type i=0; i<bound; ++i) {
            optr[i] = op(lunit * rptr[i]);
        }
    }

    template <typename C, typename Fn>
    LA_INLINE_ALWAYS static void
    linear_fwd_mul_right(pointer<C> optr, const_pointer<C> rptr, const_reference<C> lunit, index_type bound, Fn op) noexcept
    {
        for (index_type i=0; i<bound; ++i) {
            optr[i] += op(lunit * rptr[i]);
        }
    }

    template<typename C, typename Fn>
    LA_INLINE_ALWAYS static void impl_0bd(pointer<C> tile,
                                          const_reference<C> lhs_unit,
                                          const_pointer<C> rhs_ptr,
                                          IDIMN stride,
                                          IDIMN ibound,
                                          IDIMN jbound,
                                          Fn op) noexcept
    {
        constexpr auto tile_width = helper_type<C>::tile_width;
        for (IDIMN i = 0; i < ibound; ++i) {
            for (IDIMN j = 0; j < jbound; ++j) {
                tile[i * tile_width + j] += op(lhs_unit * rhs_ptr[i * stride + j]);
            }
        }
    }

    template<typename C, typename Fn>
    LA_INLINE_ALWAYS static void impl_db0(pointer<C> tile,
                                          const_pointer<C> lhs_ptr,
                                          const_reference<C> rhs_unit,
                                          IDIMN stride,
                                          IDIMN ibound,
                                          IDIMN jbound,
                                          Fn op) noexcept
    {
        constexpr auto tile_width = helper_type<C>::tile_width;
        for (IDIMN i = 0; i < ibound; ++i) {
            for (IDIMN j = 0; j < jbound; ++j) {
                tile[i * tile_width + j] += op(lhs_ptr[i * stride + j] * rhs_unit);
            }
        }
    }

    template<typename C, typename Fn>
    LA_INLINE_ALWAYS static void impl_mid(pointer<C> tile,
                                          const_pointer<C> lhs_tile,
                                          const_pointer<C> rhs_tile,
                                          const int* perm,
                                          Fn op) noexcept
    {
        constexpr auto tile_width = helper_type<C>::tile_width;
        for (IDIMN i = 0; i < tile_width; ++i) {
            auto pi = perm[i];
            for (IDIMN j = 0; j < tile_width; ++j) {
                tile[i * tile_width + j] += op(lhs_tile[pi] * rhs_tile[j]);
            }
        }
    }

    template<typename C, typename Fn>
    LA_INLINE_ALWAYS static void impl_lb1(pointer<C> tile,
                                          const_reference<C> lhs_val,
                                          const_pointer<C> rhs_tile,
                                          Fn op,
                                          IDIMN i) noexcept
    {
        constexpr auto tile_width = helper_type<C>::tile_width;
        for (IDIMN j = 0; j < tile_width; ++j) {
            tile[i * tile_width + j] += op(lhs_val * rhs_tile[j]);
        }
    }

    template<typename C, typename Fn>
    LA_INLINE_ALWAYS static void impl_1br(pointer<C> tile,
                                          const_pointer<C> lhs_tile,
                                          const_reference<C> rhs_val,
                                          const int* perm,
                                          Fn op,
                                          IDIMN j) noexcept
    {
        constexpr auto tile_width = helper_type<C>::tile_width;
        for (IDIMN i = 0; i < tile_width; ++i) {
            tile[i * tile_width + j] += op(lhs_tile[perm[i]] * rhs_val);
        }
    }

    template<typename C, typename Fn>
    LA_INLINE_ALWAYS static void impl_ulmd(pointer<C> tile,
                                           const_pointer<C> lhs_fwd_ptr,
                                           const_reference<C> rhs_val,
                                           Fn op,
                                           IDIMN j,
                                           IDIMN ibound,
                                           IDIMN stride) noexcept
    {
        constexpr auto tile_width = helper_type<C>::tile_width;
        for (IDIMN i = 0; i < ibound; ++i) {
            tile[i * tile_width + j] += op(lhs_fwd_ptr[i * stride] * rhs_val);
        }
    }

    template<typename C, typename Fn>
    LA_INLINE_ALWAYS static void impl_ulmd_1l(pointer<C> tile,
                                              const_pointer<C> lhs_fwd_ptr,
                                              const_pointer<C> rhs_fwd_ptr,
                                              Fn op,
                                              IDIMN lhs_stride,
                                              IDIMN ibound,
                                              IDIMN jbound) noexcept
    {
        constexpr auto tile_width = helper_type<C>::tile_width;
        for (IDIMN i = 0; i < ibound; ++i) {
            for (IDIMN j = 0; j < jbound; ++j) {
                tile[i * tile_width + j] += op(lhs_fwd_ptr[i * lhs_stride] * rhs_fwd_ptr[j]);
            }
        }
    }

protected:

    template <typename C, typename Fn>
    static void impl_outer_cases_zip(helper_type<C>& helper, Fn op) noexcept
    {
        constexpr auto tile_width = helper_type<C>::tile_width;
        const auto out_deg = helper.out_degree();
        const_pointer<C> lptr = helper.left_read_tile_ptr();
        const_pointer<C> rptr = helper.right_read_tile_ptr();
        pointer<C> optr = helper.fwd_write(out_deg);
        const_pointer<C> rsrc = helper.right_fwd_read(out_deg);

        const_reference<C> lunit = helper.left_unit();
        const_reference<C> runit = helper.right_unit();

        const auto bound = helper.num_subtiles*tsi::powers[out_deg-helper.tile_letters];

        if (helper.is_inplace()) {
            for (index_type i=0; i<bound; ++i, rsrc+=tile_width, optr+=tile_width) {
                helper.read_right_tile(rsrc);
                helper.read_left_tile(optr);
                linear_fwd_zipped_mul_inplace<C>(optr, lptr, rptr, lunit, runit, tile_width, op);
            }
        } else {
            const_pointer<C> lsrc = helper.left_fwd_read(out_deg);

            for (index_type i=0; i<bound; ++i, lsrc+=tile_width, rsrc+=tile_width, optr+=tile_width) {
                helper.read_left_tile(lsrc);
                helper.read_right_tile(rsrc);
                linear_fwd_zipped_mul<C>(optr, lptr, rptr, lunit, runit, tile_width, op);
            }
        }
    }

    template <typename C, typename Fn>
    static void impl_outer_left_only(helper_type<C>& helper, Fn op)
    {
        constexpr auto tile_width = helper_type<C>::tile_width;
        const auto out_deg = helper.out_degree();
        const_pointer<C> lptr = helper.left_read_tile_ptr();
        const_reference<C> runit = helper.right_unit();
        pointer<C> optr = helper.fwd_write(out_deg);

        const auto bound = helper.num_subtiles * tsi::powers[out_deg - helper.tile_letters];

        if (helper.is_inplace()) {
            for (index_type i=0; i<bound; ++i, optr+=tile_width) {
                helper.read_left_tile(optr);
                linear_fwd_mul_inplace_left<C>(optr, lptr, runit, tile_width, op);
            }
        } else {
            const_pointer<C> lsrc = helper.left_fwd_read(out_deg);

            for (index_type i=0; i<bound; ++i, lsrc+=tile_width, optr+=tile_width) {
                helper.read_left_tile(lsrc);
                linear_fwd_mul_left<C>(optr, lptr, runit, tile_width, op);
            }
        }
    }

    template <typename C, typename Fn>
    static void impl_outer_right_only(helper_type<C>& helper, Fn op)
    {
        constexpr auto tile_width = helper_type<C>::tile_width;
        const auto out_deg = helper.out_degree();
        const_pointer<C> rptr = helper.right_read_tile_ptr();
        const_reference<C> lunit = helper.left_unit();
        pointer<C> optr = helper.fwd_write(out_deg);

        const auto bound = helper.num_subtiles * tsi::powers[out_deg - helper.tile_letters];

        if (helper.is_inplace()) {
            for (index_type i=0; i<bound; ++i, optr+=tile_width) {
                helper.read_right_tile(optr);
                linear_fwd_mul_inplace_left<C>(optr, rptr, lunit, tile_width, op);
            }
        } else {
            const_pointer<C> rsrc = helper.right_fwd_read(out_deg);

            for (index_type i=0; i<bound; ++i, rsrc+=tile_width, optr+=tile_width) {
                helper.read_right_tile(rsrc);
                linear_fwd_mul_left<C>(optr, rptr, lunit, tile_width, op);
            }
        }
    }

    template <typename C, typename Fn>
    static void impl_outer_cases(helper_type<C>& helper, Fn op) noexcept
    {
        bool left_valid = helper.lhs_degree() >= helper.out_degree()
                && helper.left_unit() != C::zero;
        bool right_valid = helper.rhs_degree() >= helper.out_degree()
                && helper.right_unit() != C::zero;

        if (left_valid && right_valid) {
            impl_outer_cases_zip(helper, op);
        } else if (left_valid) {
            impl_outer_left_only(helper, op);
        } else {
            impl_outer_right_only(helper, op);
        }
    }

    template<typename Coeffs, typename Fn>
    static void
    impl_lhs_small(helper_type<Coeffs>& helper,
                   Fn op,
                   IDEG out_deg,
                   IDEG lhs_deg,
                   IDIMN k,
                   IDIMN subtile_i,
                   IDIMN subtile_j) noexcept
    {
        const auto rhs_deg = out_deg - lhs_deg;
        constexpr auto tile_width = helper_type<Coeffs>::tile_width;
        constexpr auto tile_letters = helper_type<Coeffs>::tile_letters;

        const auto ibound = helper.boundary_subtile(subtile_i) ? Width % tile_width : tile_width;

        assert(1 <= lhs_deg && lhs_deg <= tile_letters);
        for (IDIMN i = 0; i < ibound; ++i) {
            const auto split = helper.split_key(tile_letters - lhs_deg, subtile_i * tile_width + i);
            const auto& left_val = *helper.left_fwd_read(lhs_deg, split.first);
            helper.read_right_tile(rhs_deg,
                                   helper.combine_keys(out_deg - 2 * tile_letters, split.second, k),
                                   subtile_j);
            impl_lb1<Coeffs>(helper.out_tile_ptr(),
                             left_val,
                             helper.right_read_tile_ptr(),
                             op,
                             i);
        }
    }

    template<typename Coeffs, typename Fn>
    static void
    impl_mid_cases_reverse(helper_type<Coeffs>& helper,
                           Fn op,
                           IDEG out_deg,
                           IDEG lhs_deg,
                           IDIMN k,
                           IDIMN subtile_i,
                           IDIMN subtile_j) noexcept
    {
        constexpr auto tile_letters = helper_type<Coeffs>::tile_letters;

        const auto rhs_deg = out_deg - lhs_deg;
        assert(tile_letters <= lhs_deg && lhs_deg <= out_deg - tile_letters);
        assert(tile_letters <= rhs_deg && rhs_deg <= out_deg - tile_letters);
        const auto lhs_split = lhs_deg - tile_letters;
        const auto rhs_split = rhs_deg - tile_letters;
        assert(lhs_split + rhs_split == out_deg - 2 * tile_letters);

        const auto split = helper.split_key(rhs_split, k);
        helper.read_left_tile(lhs_deg, helper.reverse_key(lhs_split, split.first), subtile_i);
        helper.read_right_tile(rhs_deg, split.second, subtile_j);
        impl_mid<Coeffs>(helper.out_tile_ptr(),
                         helper.left_read_tile_ptr(),
                         helper.right_read_tile_ptr(),
                         helper.reverser(),
                         op);
    }

    template<typename Coeffs, typename Fn>
    static void
    impl_mid_cases_no_reverse(helper_type<Coeffs>& helper,
                              Fn op,
                              IDEG out_deg,
                              IDEG lhs_deg,
                              IDIMN k,
                              IDIMN subtile_i,
                              IDIMN subtile_j) noexcept
    {
        constexpr auto tile_width = helper_type<Coeffs>::tile_width;
        constexpr auto tile_letters = helper_type<Coeffs>::tile_letters;

        const auto rhs_deg = out_deg - lhs_deg;
        const auto lhs_split = lhs_deg - tile_letters;
        const auto rhs_split = rhs_deg - tile_letters;
        assert(lhs_split + rhs_split == out_deg - 2 * tile_letters);

        const auto split = helper.split_key(rhs_split, k);
        const auto lhs_stride = tsi::powers[lhs_split];
        impl_ulmd_1l<Coeffs>(helper.out_tile_ptr(),
                             helper.left_fwd_read(lhs_deg, split.first + subtile_i * tile_width * lhs_stride),
                             helper.right_fwd_read_ptr(rhs_deg, split.second, subtile_j),
                             op,
                             lhs_stride,
                             helper.boundary_subtile(subtile_i) ? Width % tile_width : tile_width,
                             helper.boundary_subtile(subtile_j) ? Width % tile_width : tile_width);
    }

    template<typename Coeffs, typename Fn>
    static void
    impl_rhs_small_reverse(helper_type<Coeffs>& helper,
                           Fn op,
                           IDEG out_deg,
                           IDEG lhs_deg,
                           IDIMN k_reverse,
                           IDIMN subtile_i,
                           IDIMN subtile_j) noexcept
    {
        constexpr auto tile_width = helper_type<Coeffs>::tile_width;
        constexpr auto tile_letters = helper_type<Coeffs>::tile_letters;

        const auto rhs_deg = out_deg - lhs_deg;
        assert(out_deg - 2 * tile_letters < lhs_deg && lhs_deg < out_deg);
        assert(1 <= rhs_deg && rhs_deg < tile_letters);
        const auto split_left_letters = tile_letters - rhs_deg;
        assert(0 < split_left_letters && split_left_letters < tile_letters);
        assert(lhs_deg == out_deg - 2 * tile_letters + split_left_letters + tile_letters);

        for (IDIMN j = 0; j < tile_width; ++j) {
            const auto split = helper.split_key(rhs_deg, subtile_j * tile_width + j);
            const auto& right_val = *helper.right_fwd_read(rhs_deg, split.second);
            helper.read_left_tile(lhs_deg,
                                  helper.combine_keys(split_left_letters, helper.reverse_key(split_left_letters, split.first), k_reverse),
                                  subtile_i);
            impl_1br<Coeffs>(helper.out_tile_ptr(),
                             helper.left_read_tile_ptr(),
                             right_val, helper.reverser(),
                             op,

                             j);
        }
    }

    template<typename Coeffs, typename Fn>
    static void
    impl_rhs_small_no_reverse(helper_type<Coeffs>& helper,
                              Fn op,
                              IDEG out_deg,
                              IDEG lhs_deg,
                              IDIMN k,
                              IDIMN subtile_i,
                              IDIMN subtile_j) noexcept
    {
        constexpr auto tile_width = helper_type<Coeffs>::tile_width;
        constexpr auto tile_letters = helper_type<Coeffs>::tile_letters;

        const auto rhs_deg = out_deg - lhs_deg;
        const auto split_left_letters = tile_letters - rhs_deg;
        const auto lhs_stride = tsi::powers[lhs_deg - tile_letters];
        for (IDIMN j = 0; j < tile_width; ++j) {
            const auto split = helper.split_key(rhs_deg, subtile_j * tile_width + j);
            const auto& right_val = *helper.right_fwd_read(rhs_deg, split.second);
            const auto lhs_key = helper.combine_keys(split_left_letters, k, split.first);

            impl_ulmd<Coeffs>(helper.out_tile_ptr(),
                              helper.left_fwd_read(lhs_deg, lhs_key + subtile_i * tile_width * lhs_stride),
                              right_val,
                              op,
                              j,
                              helper.boundary_subtile(subtile_i) ? Width % tile_width : tile_width,
                              lhs_stride);
        }
    }

    template<typename Coeffs, typename Fn>
    void impl_common(helper_type<Coeffs>& helper, Fn op) const
    {
        constexpr auto tile_width = helper_type<Coeffs>::tile_width;
        constexpr auto tile_letters = helper_type<Coeffs>::tile_letters;
        constexpr auto num_subtiles = helper_type<Coeffs>::num_subtiles;

        const auto max_degree = helper.out_degree();
        const auto old_lhs_deg = helper.lhs_degree();
        const auto rhs_max_deg = helper.rhs_degree();

        std::array<const_pointer<Coeffs>, Depth + 1> left_reads = {};
        std::array<const_pointer<Coeffs>, Depth + 1> right_reads = {};

        left_reads[0] = helper.left_fwd_read(0, 0);
        right_reads[0] = helper.right_fwd_read(0, 0);

        impl_outer_cases(helper, op);

        for (IDEG out_deg = max_degree; out_deg > 2 * tile_letters; --out_deg) {
            const auto mid_deg = out_deg - 2 * tile_letters;
            const auto mid_end = out_deg - tile_letters;
            const auto stride = static_cast<IDIMN>(tsi::powers[out_deg - tile_letters]);

            auto lhs_deg_min = std::max(IDEG(1), out_deg - rhs_max_deg);
            auto lhs_deg_max = std::min(out_deg - 1, old_lhs_deg);

            for (IDIMN k = 0; k < static_cast<IDIMN>(tsi::powers[mid_deg]); ++k) {
                auto k_reverse = helper.reverse_key(mid_deg, k);

                for (IDIMN subtile_i = 0; subtile_i < num_subtiles; ++subtile_i) {
                    const auto ibound = helper.boundary_subtile(subtile_i) ? Width % tile_width : tile_width;
                    for (IDIMN subtile_j = 0; subtile_j < num_subtiles; ++subtile_j) {
                        const auto jbound = helper.boundary_subtile(subtile_j) ? Width % tile_width : tile_width;

                        helper.reset_tile(out_deg, k, k_reverse, subtile_i, subtile_j);

                        // Setup read pointers

                        for (IDEG i = std::max(tile_letters, lhs_deg_min); i <= std::min(out_deg - tile_letters, lhs_deg_max); ++i) {
                            auto split = helper.split_key(out_deg - i - tile_letters, k);
                            auto lkey = helper.reverse_key(i - tile_letters, split.first);
                            auto rkey = split.second;
                            left_reads[i] = helper.left_rev_read_ptr(i, lkey, subtile_i);
                            right_reads[out_deg - i] = helper.right_fwd_read_ptr(out_deg - i, rkey, subtile_j);
                        }
                        if (out_deg <= old_lhs_deg) {
                            left_reads[out_deg] = helper.left_fwd_read_ptr(out_deg, k, subtile_i);
                        }
                        if (out_deg <= rhs_max_deg) {
                            right_reads[out_deg] = helper.right_fwd_read_ptr(out_deg, k, subtile_j);
                        }

                        if (out_deg < max_degree) {
                            const auto& rhs_unit = helper.right_unit();
                            if (out_deg <= old_lhs_deg && rhs_unit != Coeffs::zero) {
                                impl_db0<Coeffs>(helper.out_tile_ptr(), left_reads[out_deg], rhs_unit, stride, ibound, jbound, op);
                            }
                            const auto& lhs_unit = helper.left_unit();
                            if (out_deg <= rhs_max_deg && lhs_unit != Coeffs::zero) {
                                impl_0bd<Coeffs>(helper.out_tile_ptr(), lhs_unit, right_reads[out_deg], stride, ibound, jbound, op);
                            }
                        }

                        for (IDEG lhs_deg = lhs_deg_min; lhs_deg <= lhs_deg_max; ++lhs_deg) {
                            if (lhs_deg < tile_letters) {
                                impl_lhs_small(helper, op, out_deg, lhs_deg, k, subtile_i, subtile_j);
                            }
                            else if (lhs_deg <= mid_end && lhs_deg < old_lhs_deg) {
                                helper.read_left_tile(left_reads[lhs_deg], helper.subtile_bound(subtile_i));
                                helper.read_right_tile(right_reads[out_deg - lhs_deg], helper.subtile_bound(subtile_j));
                                impl_mid<Coeffs>(helper.out_tile_ptr(), helper.left_read_tile_ptr(), helper.right_read_tile_ptr(), helper.reverser(), op);
                                //                                impl_mid_cases_reverse(helper, op, out_deg, lhs_deg, k, subtile_i, subtile_j);
                            }
                            else if (lhs_deg <= mid_end && lhs_deg == old_lhs_deg) {
                                impl_mid_cases_no_reverse(helper, op, out_deg, lhs_deg, k, subtile_i, subtile_j);
                            }
                            else if (lhs_deg > mid_end && lhs_deg < old_lhs_deg) {
                                impl_rhs_small_reverse(helper, op, out_deg, lhs_deg, k_reverse, subtile_i, subtile_j);
                            }
                            else if (lhs_deg > mid_end) {
                                impl_rhs_small_no_reverse(helper, op, out_deg, lhs_deg, k, subtile_i, subtile_j);
                            }
                            else {
                                BOOST_UNREACHABLE_RETURN()
                            }
                        }


                        helper.write_tile(out_deg, k, k_reverse, subtile_i, subtile_j);
                    }// subtile_j
                }    // subtile_i
            }        // k
        }            // out_deg
    }

    template<typename Coeffs, typename Fn>
    void fma_impl(helper_type<Coeffs>& helper, Fn op, DEG max_degree) const
    {
        impl_common(helper, op);
        base::fma_impl(helper, op, std::min(max_degree, DEG(2 * helper_type<Coeffs>::tile_letters)));
    }

    template<typename Coeffs, typename Fn>
    void multiply_inplace_impl(helper_type<Coeffs>& helper, Fn op, DEG max_degree) const
    {
        impl_common(helper, op);
        base::multiply_inplace_impl(helper, op, std::min(IDEG(max_degree), 2 * helper_type<Coeffs>::tile_letters));
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
        if (max_degree <= 2 * helper_type<Coeffs>::tile_letters) {
            base::fma(out, lhs, rhs, op, max_degree, orig);
            return;
        }

        if (lhs.dimension() != 0 && rhs.dimension() != 0) {
            helper_type<Coeffs> helper(out, lhs, rhs, max_degree);
            fma_impl(helper, op, helper.out_degree());
            if (helper.out_degree() > 0) {
                auto update_deg = std::min(helper.out_degree() - 1, 2 * helper_type<Coeffs>::tile_letters);
                base::update_reverse_data(out, update_deg);
            }
        }
    }

    template<typename Basis, typename Coeffs, typename Fn, typename OriginalVectors>
    void multiply_inplace(vectors::dense_vector<Basis, Coeffs>& lhs,
                          const vectors::dense_vector<Basis, Coeffs>& rhs,
                          Fn op, DEG max_degree, OriginalVectors& orig) const
    {
        //        std::cout << "BEFORE " << lhs << '\n';
        if (max_degree <= 2 * helper_type<Coeffs>::tile_letters) {
            base::multiply_inplace(lhs, rhs, op, max_degree, orig);
            return;
        }

        if (rhs.dimension() != 0) {
            helper_type<Coeffs> helper(lhs, rhs, max_degree);
            multiply_inplace_impl(helper, op, helper.out_degree());
            if (helper.out_degree() > 0) {
                auto update_deg = std::min(helper.out_degree() - 1, 2 * helper_type<Coeffs>::tile_letters);
                base::update_reverse_data(lhs, update_deg);
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
using free_tensor_multiplication = tiled_free_tensor_multiplication<Width, Depth>;

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
class free_tensor
    : public algebra<
              free_tensor_basis<n_letters, max_degree>,
              Coeff,
              free_tensor_multiplication<n_letters, max_degree>,
              VectorType,
              free_tensor<Coeff, n_letters, max_degree, VectorType, Args...>,
              Args...>
{
    typedef free_tensor_multiplication<n_letters, max_degree> multiplication_t;

    using base = algebra<
            free_tensor_basis<n_letters, max_degree>,
            Coeff,
            free_tensor_multiplication<n_letters, max_degree>,
            VectorType,
            free_tensor,
            Args...>;

    template<template<typename, typename, typename...> class VT, typename... VArgs>
    static void resize_for_degree(free_tensor<Coeff, n_letters, max_degree, VT, VArgs...>& arg, DEG degree)
    {}

    template<typename... VArgs>
    static void resize_for_degree(free_tensor<Coeff, n_letters, max_degree, ::alg::vectors::dense_vector, VArgs...>& arg, DEG degree)
    {
        arg.base_vector().resize_to_degree(degree);
    }

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

    using basis_type = free_tensor_basis<n_letters, max_degree>;
    using key_type = typename basis_type::KEY;
    using scalar_type = typename Coeff::S;

public:
    using base::base;

    explicit free_tensor(typename boost::call_traits<scalar_type>::param_type s)
        : base(key_type{}, s)
    {}

    template<typename Letter, typename Scalar>
    explicit free_tensor(Letter let, Scalar sca)
        : base(base::basis.keyofletter(LET(let)), scalar_type(sca))
    {}

    /// Computes the truncated exponential of a free_tensor instance.
    friend free_tensor exp(const free_tensor& arg)
    {
        // Computes the truncated exponential of arg
        // 1 + arg + arg^2/2! + ... + arg^n/n! where n = max_degree
        KEY kunit;
        free_tensor result(kunit);

        resize_for_degree(result, max_degree);

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

        resize_for_degree(result, max_degree);

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

        resize_for_degree(*this, max_degree);

        if (unit_elt != x.end() && unit_elt->value() != typename free_tensor::SCALAR(0)) {
            x.erase(unit_elt);
        }

        for (DEG i = max_degree; i >= 1; --i) {
            this->mul_scal_div(x, typename free_tensor::SCALAR(i), max_degree - i + 1);
            *this += original;
        }

        return *this;
    }
    /// Computes the antipode of a free_tensor instance.
    inline friend free_tensor antipode(const free_tensor& arg)
    {
        // Get the trait to access the storage tag, although it now occurs to me that we already had
        // the vector type as a template argument, so we might be able to dispatch off that instead
        // of a tag. But the tag will do for now.

        // Now use tagged dispatch to pick the correct implementation
        free_tensor result;
        dtl::tiled_inverse_operator<n_letters, max_degree, Coeff, dtl::default_signer>::apply(arg.base_vector(), result.base_vector());

        return result;
    }

    /// Computes the truncated logarithm of a free_tensor instance.
    friend free_tensor log(const free_tensor& arg)
    {
        // Computes the truncated log of arg up to degree max_degree
        // The coef. of the constant term (empty word in the monoid) of arg
        // is forced to 1.
        // log(arg) = log(1+x) = x - x^2/2 + ... + (-1)^(n+1) x^n/n.
        // max_degree must be > 0

        KEY kunit;
        free_tensor tunit(kunit);
        free_tensor x(arg);
        auto it = x.find(kunit);
        if (it != x.end()) {
            x.erase(it);
        }
        free_tensor result;

        for (DEG i = max_degree; i >= 1; --i) {
            if (i % 2 == 0) {
                result.sub_scal_div(tunit, typename Coeff::Q(i));
            }
            else {
                result.add_scal_div(tunit, typename Coeff::Q(i));
            }
            result *= x;
        }

        return result;
    }

    /// Computes the truncated inverse of a free_tensor instance.
    friend free_tensor
    inverse(const free_tensor& arg)
    {
        // Computes the truncated inverse of arg up to degree max_degree
        // An exception is thrown if the leading term is zero.
        // the module assumes
        // (a+x)^(-1) = (a(1+x/a))^(-1)
        //  = a^(-1)(1 - x/a + x^2/a^2 + ... + (-1)^(n) x^n/a^n)
        // = a^(-1) - x/a*[a^(-1)(1 - x/a + x^2/a^2 + ... + (-1)^(n)
        // x^(n-1)/a^(n-1)))]. S_n = a^(-1) + z S_{n-1}; z = - x/a ; S_0 = a^(-1)
        // max_degree must be > 0

        KEY kunit;
        scalar_type a(0);
        free_tensor x, z(a);

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
        free_tensor free_tensor_a_inverse(scalar_type(1) / a), result(free_tensor_a_inverse);
        resize_for_degree(result, max_degree);

        // z := - x/a
        z.sub_scal_div(x, a);
        // the iteration
        for (DEG i = 0; i != max_degree; ++i) {
            auto tmp = z * result;
            result = free_tensor_a_inverse + z * result;
        }
        return result;
    }

    /// Computes the truncated inverse of a free_tensor instance.
    friend free_tensor
    reflect(const free_tensor& arg)
    {
        return antipode(arg);
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
        {
            return m_valid;
        }

        constexpr operator bool() const noexcept
        {
            return valid() && storage_type::size() > 0;
        }

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

public:
    void construct_reverse_data(DEG degree)
    {
        if (degree > 0) {
            if (m_reverse_data.size() < tsi::degree_sizes[degree]) {
                m_reverse_data.resize(tsi::degree_sizes[degree]);
            }
            assert(m_reverse_data.size() >= tsi::degree_sizes[degree]);
            ::alg::dtl::tiled_inverse_operator<Width, Depth, Coeffs, ::alg::dtl::non_signing_signer> t;
            t(m_data.begin(), m_reverse_data.begin(), degree);
            assert(*m_data.begin() == *m_reverse_data.begin());
        }
        m_reverse_data.validate();
    }

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
        resize_to_dimension(idx + 1);
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
        if (key.size() < m_degree) {
            m_reverse_data.invalidate();
        }
        return value(idx);
    }

private:
    void do_reserve(::alg::basis::dtl::resize_info info)
    {
        m_data.reserve(info.size);
        m_dimension = info.dimension;
        m_degree = info.degree;
        if (m_reverse_data && m_degree > 0) {
            m_reverse_data.reserve(BASIS::start_of_degree(m_degree - 1));
        }
    }

    void do_resize(::alg::basis::dtl::resize_info info)
    {
        m_data.resize(info.size);
        m_dimension = info.dimension;
        m_degree = info.degree;
        if (m_reverse_data && m_degree > 0) {
            m_reverse_data.resize(BASIS::start_of_degree(m_degree - 1));
        }
    }

public:
    /// Reserve to dimension
    void reserve_to_dimension(const DIMN dim)
    {
        if (dim > m_dimension) {
            auto info = basis_traits::next_resize_dimension(base_vector_type::basis, dim, m_degree);
            do_reserve(info);
        }
        assert(m_data.size() == m_dimension);
    }

    /// Reserve to degree
    void reserve_to_degree(const DEG deg)
    {
        if (deg > m_degree) {
            auto info = basis_traits::next_resize_dimension(base_vector_type::basis, 0, deg);
            do_reserve(info);
        }
        assert(m_data.size() == m_dimension);
    }

    DIMN resize_for_key(key_type key)
    {
        auto info = basis_traits::key_resize_dimension(base_vector_type::basis, key);
        if (info.dimension > m_dimension) {
            do_resize(info);
        }
        assert(m_data.size() == m_dimension);
        return basis_traits::key_to_index(base_vector_type::basis, key);
    }

    /// Resize to dimension
    void resize_to_dimension(DIMN dim)
    {
        auto info = basis_traits::next_resize_dimension(base_vector_type::basis, dim, m_degree);
        if (info.dimension > m_dimension) {
            do_resize(info);
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
            do_resize(info);
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

    const reverse_storage& reverse_data() const noexcept
    {
        return m_reverse_data;
    }
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
        if (m_reverse_data && key.size() < m_degree) {
            const auto idx = basis.key_to_index(key.reverse());
            m_reverse_data[idx] = Coeffs::zero;
        }
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
        operator[](key) += scalar_type(s);
        if (m_reverse_data && key.size() < m_degree) {
            const auto idx = basis.key_to_index(key.reverse());
            m_reverse_data[idx] += scalar_type(s);
        }
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
            for (DIMN i = 0; i < arg.m_reverse_data.size(); ++i) {
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
        for (DIMN i = 0; i < arg.m_reverse_data.size(); ++i) {
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

        DIMN i = 0;
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

            for (; i < rminmax.first; ++i) {
                result.m_reverse_data.emplace(i, op(lhs.m_reverse_data[i], rhs.m_reverse_data[i]));
            }
            for (; i < lhs_rsize; ++i) {
                result.m_reverse_data.emplace(i, op(lhs.m_reverse_data[i], Coeffs::zero));
            }
            for (; i < rhs_rsize; ++i) {
                result.m_reverse_data.emplace(i, op(Coeffs::zero, rhs.m_reverse_data[i]));
            }
            result.m_reverse_data.validate();
        }
        else {
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

        if (minmax.second > lhs.m_dimension) {
            lhs.reserve_to_dimension(minmax.second);
        }
        const auto size_fwd = minmax.first;

        DIMN i = 0;

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

            i = 0;
            for (; i < rminmax.first; ++i) {
                lhs.m_reverse_data[i] = op(lhs.m_reverse_data[i],
                                           rhs.m_reverse_data[i]);
            }
            for (; i < lhs_rsize; ++i) {
                lhs.m_reverse_data[i] = op(lhs.m_reverse_data[i], Coeffs::zero);
            }
            for (; i < rhs_rsize; ++i) {
                lhs.m_reverse_data.emplace(i, op(Coeffs::zero, rhs.m_reverse_data[i]));
            }
            lhs.m_reverse_data.validate();
        }
        else {
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

        for (auto i = mid_eq.first; i < m_dimension; ++i) {
            if (m_data[i] != Coeffs::zero) {
                return false;
            }
        }

        for (auto i = mid_eq.first; i < rhs.m_dimension; ++i) {
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
