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

#include <algorithm>
#include <unordered_set>

#include <boost/align/aligned_alloc.hpp>
#include <boost/align/aligned_allocator.hpp>
#include <boost/align/assume_aligned.hpp>
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

#include "detail/unpacked_tensor_word.h"
#include "detail/platform.h"
#include "half_shuffle_tensor_basis.h"
#include "tensor_basis.h"


#define LA_RESTRICT __restrict

#if BOOST_COMP_GNUC || BOOST_COMP_CLANG
#define LA_INLINE_ALWAYS __attribute__((always_inline))
#else
#define LA_INLINE_ALWAYS
#endif


#define LA_ALIGNAS(BYTES) alignas(BYTES)
#ifndef LA_CACHELINE_BYTES
#define LA_CACHELINE_BYTES 64
#endif

#ifndef LIBALGEBRA_L1_CACHE_SIZE
#define LIBALGEBRA_L1_CACHE_SIZE 32768// 32KB should be fairly standard
#endif

#define LA_DEFAULT_TILE_PARAM(WIDTH, DEPTH, SCALAR) ::alg::dtl::tensor_tile_letters_helper<WIDTH, DEPTH, LIBALGEBRA_L1_CACHE_SIZE / (2 * sizeof(SCALAR))>::num_letters


// Define LIBALGEBRA_TM_UPDATE_REVERSE_INLINE to allow writing of reverse data during multiplication
#define LIBALGEBRA_TM_UPDATE_REVERSE_INLINE


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

template<DEG Width, DEG Depth, DIMN TargetSize, bool WidthSquareFits = (Width * Width < TargetSize)>
struct tensor_tile_letters_helper {
    static constexpr IDEG log_target = IDEG(::alg::integer_maths::logN(TargetSize, Width));
    static constexpr IDEG num_letters_unlimited = std::min(IDEG(Depth) / 2, (log_target >= 2) ? log_target / 2 : -1);
#ifdef LIBALGEBRA_MAX_TILE_LETTERS
    static constexpr IDEG num_letters = std::min(
            IDEG(LIBALGEBRA_MAX_TILE_LETTERS),
            num_letters_unlimited);
#else
    static constexpr IDEG num_letters = num_letters_unlimited;
#endif
};

template <DEG Depth, DIMN T, bool B>
struct tensor_tile_letters_helper<1, Depth, T, B> {
    static constexpr IDEG num_letters = 0;
};


template<DEG Width, DEG Depth, DIMN TargetSize>
struct tensor_tile_letters_helper<Width, Depth, TargetSize, false> {
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
//        data = new S[Size] {};
        std::uninitialized_fill_n(data, Size, S());
    }

    ~data_tile()
    {
        boost::alignment::aligned_free(data);
//        delete[] data;
    }
};

template<DEG Width, DEG Depth, IDEG TileLetters = 0>
struct tile_details {
    static constexpr IDEG tile_letters = (Depth > 1) ? (TileLetters > 0) ? TileLetters : 1 : 0;
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
    void read_tile_impl(const_pointer LA_RESTRICT in_p, index_type lhs_stride, index_type ibound, index_type jbound) noexcept
    {
        for (index_type i = 0; i < ibound; ++i) {
            for (index_type j = 0; j < jbound; ++j) {
                tile.data[i * tile_width + j] = in_p[i * lhs_stride + j];
            }
        }
    }

    template <index_type IBound, index_type JBound, index_type OutStride=JBound>
    inline LA_INLINE_ALWAYS void read_tile_impl(const_pointer LA_RESTRICT src, index_type in_stride)
    {
        for (index_type i = 0; i < IBound; ++i) {
            for (index_type j = 0; j < IBound; ++j) {
                tile.data[i * OutStride + j] = src[i * in_stride + j];
            }
        }
    }


    template <index_type IBound, index_type JBound, index_type InStride=JBound>
    LA_INLINE_ALWAYS static void write_tile_assign(pointer __restrict optr, const_pointer __restrict tptr, index_type out_stride)
    {
        BOOST_ALIGN_ASSUME_ALIGNED(tptr, LA_CACHELINE_BYTES);
        for (index_type i = 0; i < IBound; ++i) {
            LA_PREFETCH_ET0(optr + (i+1)*out_stride);
            for (index_type j = 0; j < JBound; ++j) {
                optr[i*out_stride + j] = tptr[i*InStride + j];
            }
        }
    }

    LA_INLINE_ALWAYS static void write_tile_assign(pointer __restrict out_p, const_pointer __restrict tptr, index_type ibound, index_type jbound, index_type out_stride, index_type in_stride) noexcept
    {
        BOOST_ALIGN_ASSUME_ALIGNED(tptr, LA_CACHELINE_BYTES);
        for (index_type i = 0; i < ibound; ++i) {
            LA_PREFETCH_ET0(out_p + (i+1)*out_stride);

            for (index_type j = 0; j < jbound; ++j) {
                out_p[i*out_stride + j] = tptr[i*in_stride + j];
            }
        }
    }

    template <index_type IBound, index_type JBound, index_type InStride=JBound>
    LA_INLINE_ALWAYS static void write_tile_acc(pointer __restrict optr, const_pointer __restrict tptr, index_type out_stride)
    {
        BOOST_ALIGN_ASSUME_ALIGNED(tptr, LA_CACHELINE_BYTES);
        for (index_type i = 0; i < IBound; ++i) {
            LA_PREFETCH_ET0(optr + (i+1)*out_stride);
            for (index_type j = 0; j < IBound; ++j) {
                optr[i*out_stride + j] += tptr[i*InStride + j];
            }
        }
    }

    LA_INLINE_ALWAYS static void write_tile_acc(pointer __restrict optr, const_pointer __restrict tptr, index_type ibound, index_type jbound, index_type out_stride, index_type in_stride)
    {
        BOOST_ALIGN_ASSUME_ALIGNED(tptr, LA_CACHELINE_BYTES);
        for (index_type i = 0; i < ibound; ++i) {
            LA_PREFETCH_ET0(optr + (i+1)*out_stride);
            for (index_type j = 0; j < jbound; ++j) {
                optr[i*out_stride + j] += tptr[i*in_stride + j];
            }
        }
    }

    void write_tile_reverse_impl(pointer out_p, index_type ibound, index_type jbound, index_type out_stride, index_type in_stride) const noexcept
    {
        using perm = reversing_permutation<Width, tile_info::tile_letters>;
        for (index_type i = 0; i < ibound; ++i) {
            for (index_type j = 0; j < jbound; ++j) {
                out_p[i * out_stride + j] = tile.data[perm::permute_idx(j) * in_stride + perm::permute_idx(i)];
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
        if (tile_info::num_subtiles == 1) {
            read_tile_impl<tile_info::tile_width, tile_info::tile_width>(in_p, stride);
        } else {
            read_tile_impl(in_p, stride, subtile_bound(subtile_i), subtile_bound(subtile_j));
        }
    }

    void write_tile(pointer dst_p,
                    degree_type degree,
                    index_type index,
                    index_type subtile_i = 0,
                    index_type subtile_j = 0) const noexcept
    {
        const auto stride = tsi::powers[degree - tile_info::tile_letters];
        pointer out_p = dst_p + pointer_offset(degree, index, subtile_i, subtile_j);
        if (tile_info::num_subtiles == 1) {
            write_tile_assign<tile_info::tile_width, tile_info::tile_width>(out_p, tile.data, stride);
        } else {
            write_tile_assign(out_p, tile.data, subtile_bound(subtile_i), subtile_bound(subtile_j), stride, tile_width);
        }
    }

    void write_tile_reverse(pointer dst_p,
                            degree_type degree,
                            index_type index,
                            index_type subtile_i = 0,
                            index_type subtile_j = 0) const noexcept
    {
        const auto stride = tsi::powers[degree - tile_info::tile_letters];
        pointer out_p = dst_p + pointer_offset(degree, index, subtile_j, subtile_i);
        write_tile_reverse_impl(out_p, subtile_bound(subtile_i), subtile_bound(subtile_j), stride, tile_width);
    }

    template <size_type IIndex, size_type JIndex>
    struct permute_info
    {
        using perm = dtl::reversing_permutation<Width, tile_info::tile_letters>;
        static constexpr size_type reverse_i = perm::permute_idx(IIndex);
        static constexpr size_type reverse_j = perm::permute_idx(JIndex);
    };

    template <size_type IIndex, size_type JIndex>
    static constexpr LA_INLINE_ALWAYS void permute_loop_impl(pointer LA_RESTRICT tptr, permute_info<IIndex, JIndex> marker)
    {
        std::swap(tptr[IIndex*tile_width+JIndex], tptr[marker.reverse_j*tile_width + marker.reverse_i]);
        permute_loop_impl(tptr, permute_info<IIndex, JIndex+1>());
    }

    template <size_type IIndex>
    static constexpr LA_INLINE_ALWAYS void permute_loop_impl(pointer LA_RESTRICT tptr, permute_info<IIndex, tile_width> marker)
    {
        permute_loop_impl(tptr, permute_info<IIndex + 1, IIndex+2>());
    }

    template <size_type JIndex>
    static constexpr LA_INLINE_ALWAYS void permute_loop_impl(pointer LA_RESTRICT tptr, permute_info<tile_width, JIndex> marker)
    {}

    LA_INLINE_ALWAYS constexpr void permute_write_tile()
    {
//        using perm = dtl::reversing_permutation<Width, tile_info::tile_letters>;
        pointer tptr = tile_ptr();
        const auto* perm = reverser();
        for (index_type i = 0; i < tile_width; ++i) {
//            auto pi = perm::permute_idx(i);
            auto pi = perm[i];
            for (index_type j = i + 1; j < tile_width; ++j) {
//                auto pj = perm::permute_idx(j);
                auto pj = perm[j];
                std::swap(tptr[i * tile_width + j], tptr[pj * tile_width + pi]);
            }
        }
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
         IDEG TileLetters = LA_DEFAULT_TILE_PARAM(Width, MaxDepth, typename Coeffs::S)>
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
        static inline LA_INLINE_ALWAYS void eval(S* LA_RESTRICT dst_ptr, const S* LA_RESTRICT src_ptr, DEG current_degree)
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

    static void LA_INLINE_ALWAYS untiled_cases(scalar_type* LA_RESTRICT dst_ptr, const scalar_type* LA_RESTRICT src_ptr, DEG current_degree) noexcept
    {
        dtl::increasing_level_walker<untiled_compute, 2 * tile_info::tile_letters>::eval(dst_ptr, src_ptr, current_degree);
    }

    static void LA_INLINE_ALWAYS sign_and_permute(scalar_type* LA_RESTRICT tile, Signer& signer) noexcept
    {
        for (DIMN i = 0; i < tile_info::tile_size; ++i) {
            tile[i] = signer(tile[i]);
        }

    }

    static void permute_level_tiled(helper_type<Coeffs>& helper, DEG out_deg) noexcept
    {
        static constexpr unsigned max_middle_letters = MaxDepth > 2*tile_info::tile_letters
                ? MaxDepth - 2 * tile_info::tile_letters : 0;
        const auto mid_deg = out_deg - 2 * tile_info::tile_letters;
        Signer signer(out_deg);

        unpacked_tensor_word<Width, max_middle_letters> word;
        word.reset(mid_deg);
        for (IDIMN middle_index = 0; middle_index < IDIMN(tsi::powers[mid_deg]); ++middle_index, ++word) {
//            const auto reverse_middle_index = helper.reverse_key(mid_deg, middle_index);
            const auto reverse_middle_index = word.to_reverse_index();

            for (IDIMN subtile_i = 0; subtile_i < tile_info::num_subtiles; ++subtile_i) {
                for (IDIMN subtile_j = 0; subtile_j < tile_info::num_subtiles; ++subtile_j) {
                    helper.read_tile(out_deg, middle_index, subtile_i, subtile_j);
                    sign_and_permute(helper.tile_ptr(), signer);
                    helper.permute_write_tile();
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
        for (auto item : src) {
            auto key = item.key();
            auto deg = key.size();
            if (deg <= max_degree) {
                Signer s(deg);
                result[key.reverse()] = s(item.value());
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

template <DEG MaxDepth, typename Coeffs, typename Signer, IDEG UNUSED>
class tiled_inverse_operator<1, MaxDepth, Coeffs, Signer, UNUSED> {
public:
    using scalar_type = typename Coeffs::S;

    template<typename Vector>
    static void apply(const Vector& src, Vector& result, DEG max_degree = MaxDepth)
    {
        for (auto item : src) {
            auto key = item.key();
            auto deg = key.size();
            if (deg <= max_degree) {
                Signer s(deg);
                result[key] = s(src[key]);
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

        auto* optr = dst.as_mut_ptr();
        const auto* iptr = src.as_ptr();
        for (DEG d=0; d<=dst.degree(); ++d) {
            Signer s(d);
            optr[d] = s(iptr[d]);
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

        for (DEG d = 0; d <= curr_degree; ++d) {
            Signer s(d);
            dst_ptr[d] = s(src_ptr[d]);
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
    const_reference left_unit() const noexcept { return lhs_data == nullptr ? *out_data : *lhs_data; }
    const_reference right_unit() const noexcept { return *rhs_data; }
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
    pointer fwd_write(IDEG d, index_type offset = 0) const noexcept
    {
        assert(d >= 0 && d <= out_deg);
        assert(offset >= 0 && offset < index_type(tsi::powers[d]));
        //        return out_data + basis_type::start_of_degree(d);
        return m_out_levels[d] + offset;
    }

    std::pair<DIMN, DIMN> range_size(IDEG lhs, IDEG rhs) const noexcept
    {
        return std::pair<DIMN, DIMN>(tsi::powers[lhs], tsi::powers[rhs]);
    }
};


template<DEG Width, DEG Depth, typename Coeffs, IDEG TileLetters, IDEG WriteCacheLetters>
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

    pointer fwd_write_cache = nullptr;
    pointer rev_write_cache = nullptr;

    dtl::data_tile<const_pointer, tile_info::tile_width> p_read_queue;

    static constexpr index_type cache_letters = tile_info::tile_letters + WriteCacheLetters;
    static constexpr index_type cache_ncols = integer_maths::power(Width, cache_letters);
    static constexpr index_type cache_nrows = integer_maths::power(Width, tile_info::tile_letters);
    static constexpr index_type cache_size = cache_nrows * cache_ncols;

    static constexpr index_type num_tiles_in_cache = tile_info::num_subtiles*tile_info::num_subtiles*integer_maths::power(Width, WriteCacheLetters);
    index_type cache_count = 0;

    void allocate_write_caches()
    {
        if (WriteCacheLetters > 0) {
            if (Depth > (WriteCacheLetters + 2 * tile_info::tile_letters)) {
                fwd_write_cache = static_cast<scalar_type*>(boost::alignment::aligned_alloc(64, cache_size*sizeof(scalar_type)));
//                fwd_write_cache = new scalar_type[cache_size] {};
                std::uninitialized_fill_n(fwd_write_cache, cache_size, Coeffs::zero);
            }
#ifdef LIBALGEBRA_TM_UPDATE_REVERSE_INLINE
            if (Depth >  (1 + WriteCacheLetters + 2*tile_info::tile_letters)) {
                rev_write_cache = static_cast<scalar_type*>(boost::alignment::aligned_alloc(64, cache_size*sizeof(scalar_type)));
                //                rev_write_cache = new scalar_type[cache_size] {};
                std::uninitialized_fill_n(rev_write_cache, cache_size, Coeffs::zero);
            }
#endif
        }
    }

    void deallocate_write_caches() {
//        delete[] fwd_write_cache;
        boost::alignment::aligned_free(fwd_write_cache);
#ifdef LIBALGEBRA_TM_UPDATE_REVERSE_INLINE
//        delete[] rev_write_cache;
        boost::alignment::aligned_free(rev_write_cache);
#endif
    }


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
        allocate_write_caches();
    }

    template<typename B1, typename B2>
    tiled_free_tensor_multiplication_helper(dense_tensor<B1>& lhs, const dense_tensor<B2>& rhs, DEG max_degree)
        : base(lhs, rhs, max_degree)
    {
        setup_reverse_readwrite(lhs);
        setup_reverse_levels();
        allocate_write_caches();
    }

    ~tiled_free_tensor_multiplication_helper() {
        deallocate_write_caches();
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
    pointer rev_write(IDEG degree, IDIMN offset) noexcept
    { return reverse_write_ptr + basis_type::start_of_degree(degree) + offset; }

    const_pointer* read_queue() const noexcept { return p_read_queue.data; }

private:
    template<IDEG Letters, bool IsEven = (Letters & 1) == 0>
    struct tile_transposer {
        static constexpr IDEG half_letters = Letters / 2;
        static constexpr index_type half_tile_width = integer_maths::power(index_type(Width), half_letters);
        using half_letters_permute = reversing_permutation<Width, half_letters>;

        static constexpr void LA_INLINE_ALWAYS permute(pointer tile)
        {
            for (index_type left = 0; left < half_tile_width; ++left) {
                for (index_type middle = 0; middle < static_cast<index_type>(Width); ++middle) {
                    for (index_type right = left + 1; right < half_tile_width; ++right) {
                        auto i = (left * Width + middle) * half_tile_width + right;
                        auto j = (half_letters_permute::permute_idx(right) * Width + middle) * half_tile_width + half_letters_permute::permute_idx(left);
                        std::swap(tile[i], tile[j]);
                    }
                }
            }
        }
    };

    template<IDEG Letters>
    struct tile_transposer<Letters, true> {
        static constexpr IDEG half_letters = Letters / 2;
        static constexpr index_type half_tile_width = integer_maths::power(index_type(Width), half_letters);
        using half_letters_permute = reversing_permutation<Width, half_letters>;

        static constexpr void LA_INLINE_ALWAYS permute(pointer tile)
        {
            for (index_type left = 0; left < half_tile_width; ++left) {
                for (index_type right = left + 1; right < half_tile_width; ++right) {
                    auto i = left * half_tile_width + right;
                    auto j = half_letters_permute::permute_idx(right) * half_tile_width + half_letters_permute::permute_idx(left);
                    std::swap(tile[i], tile[j]);
                }
            }
        }
    };

    template<bool IsEven>
    struct tile_transposer<1, IsEven> {

        static constexpr void permute(pointer tile)
        {
        }
    };

public:
    LA_INLINE_ALWAYS
    constexpr void permute_left_tile() noexcept
    {
        //        using perm = dtl::reversing_permutation<Width, Letters>;
//        auto* tptr = left_read_tile.data;
//        const auto* perm = base_helper::reverser();
//        for (index_type i = 0; i < tile_width; ++i) {
//                        auto pi = perm::permute_idx(i);
//            auto pi = perm[i];
//            if (pi < i) {
//                std::swap(tptr[i], tptr[pi]);
//            }
//        }
        tile_transposer<tile_letters>::permute(left_read_tile.data);
    }

private:


    template <index_type Count>
    LA_INLINE_ALWAYS static inline void read_tile(pointer __restrict tptr, const_pointer __restrict src)
    {
//        static constexpr index_type block_size = LA_CACHELINE_BYTES / sizeof(scalar_type);
        BOOST_ALIGN_ASSUME_ALIGNED(tptr, LA_CACHELINE_BYTES);
//
//        for (index_type i=0; i<Count / block_size; ++i)
//        {
//            LA_PREFETCH_T1(src + block_size);
//            for (index_type ii=0; ii<block_size; ++ii) {
//                tptr[i*block_size + ii] = src[i*block_size + ii];
//            }
//        }

        for (index_type i=0; i < Count; ++i) {
            tptr[i] = src[i];

        }

//        if (Count < tile_width) {
            for (index_type i=Count; i<tile_width; ++i) {
                tptr[i] = Coeffs::zero;
            }
//        }
    }


    LA_INLINE_ALWAYS
    static inline void read_tile(pointer __restrict tptr, const_pointer __restrict src, index_type count)
    {
//        for (index_type i=0; i<count; ++i) {
//            tptr[i] = src[i];
//        }
        BOOST_ALIGN_ASSUME_ALIGNED(tptr, LA_CACHELINE_BYTES);

        for (index_type i=0; i<count; ++i) {
            tptr[i] = src[i];
        }




//        std::copy(src, src + count, tptr);
//        for (index_type i=count; i <tile_width; ++i) {
//            tptr[i] = Coeffs::zero;
//        }
        std::fill_n(tptr + count, tile_width-count, scalar_type(0));
    }

public:

    template <index_type Count>
    LA_INLINE_ALWAYS void read_left_tile(const_pointer src)
    {
        pointer tptr = left_read_tile.data;
        read_tile<Count>(tptr, src);
    }

    LA_INLINE_ALWAYS
    void read_left_tile(const_pointer src, IDIMN count = tile_width) noexcept
    {
        pointer tptr = left_read_tile.data;
//        if (tile_info::num_subtiles == 1 && count == tile_width) {
//            read_tile<tile_width>(tptr, src);
//        } else {
            read_tile(tptr, src, count);
//        }
    }

    template <index_type Count>
    LA_INLINE_ALWAYS void read_right_tile(const_pointer src)
    {
        pointer tptr = right_read_tile.data;
        read_tile<Count>(tptr, src);
    }

    LA_INLINE_ALWAYS
    void read_right_tile(const_pointer src, IDIMN count = tile_width) noexcept
    {
        pointer tptr = right_read_tile.data;
//        if (tile_info::num_subtiles == 1 && count == tile_width) {
//            read_tile<tile_width>(tptr, src);
//        } else {
            read_tile(tptr, src, count);
//        }
    }

    template<index_type ReadCount, index_type ElementCount, typename Ptr>
    LA_INLINE_ALWAYS void read_right_fragmented(Ptr* reads)
    {
        static_assert(ReadCount*ElementCount <= tile_width,
                      "total elements cannot exceed tile width");
        pointer tptr = right_read_tile.data;
        for (index_type i=0; i<ReadCount; ++i) {
            read_tile<ElementCount>(tptr, reads[i]);
            tptr += ElementCount;
        }
    }

    template<index_type ReadCount, index_type ElementCount, typename Ptr>
    LA_INLINE_ALWAYS void read_left_fragmented(Ptr* reads)
    {
        static_assert(ReadCount*ElementCount <= tile_width,
                      "total elements cannot exceed tile width");
        pointer tptr = left_read_tile.data;
        for (index_type i=0; i<ReadCount; ++i) {
            read_tile<ElementCount>(tptr, reads[i]);
            tptr += ElementCount;
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




private:

    static void inline LA_INLINE_ALWAYS write_small_tile(pointer out, const_pointer in, index_type ibound, index_type jbound, index_type out_stride)
    {

    }

    template <unsigned D>
    void flush_cache(const unpacked_tensor_word<Width, D>& word) {
        assert(WriteCacheLetters < word.degree());
        auto degree = word.degree();  // mid degree
        auto out_degree = degree + 2*tile_info::tile_letters;
        auto stride = static_cast<index_type>(tsi::powers[out_degree - tile_letters]);
        auto index = word.split_left_index(degree - WriteCacheLetters);

        assert(index < static_cast<index_type>(tsi::powers[degree - WriteCacheLetters]));
        assert(stride > cache_ncols);

        pointer fwd_write = base::fwd_write(out_degree, index*cache_ncols);
        LA_PREFETCH_ET0(fwd_write);
        if (base::is_inplace()) {
            base_helper::template write_tile_assign<cache_nrows, cache_ncols>(fwd_write, fwd_write_cache, stride);

#ifdef LIBALGEBRA_TM_UPDATE_REVERSE_INLINE
            if (reverse_write_ptr != nullptr && IDEG(out_degree) < base::out_deg) {
                pointer reverse_write = m_out_rev_levels[out_degree];
                reverse_write += word.split_left_reverse_index(degree - WriteCacheLetters)*cache_ncols;
                LA_PREFETCH_ET0(reverse_write);
//                base_helper::permute_write_tile();
                base_helper::template write_tile_assign<cache_nrows, cache_ncols>(reverse_write, rev_write_cache, stride);
            }
#endif
        } else {
            base_helper::template write_tile_acc<cache_nrows, cache_ncols>(fwd_write, fwd_write_cache, stride);

#ifdef LIBALGEBRA_TM_UPDATE_REVERSE_INLINE
            if (reverse_write_ptr != nullptr && IDEG(out_degree) < base::out_deg) {
                pointer reverse_write = m_out_rev_levels[out_degree];
                reverse_write += word.split_left_reverse_index(degree - WriteCacheLetters)*cache_ncols;
                LA_PREFETCH_ET0(reverse_write);
//                base_helper::permute_write_tile();
                base_helper::template write_tile_acc<cache_nrows, cache_ncols>(reverse_write, rev_write_cache, stride);
            }
#endif
        }
    }

    template <unsigned D>
    void advance_counter(const unpacked_tensor_word<Width, D>& word) {
        ++cache_count;
        if (cache_count == num_tiles_in_cache) {
            /*
             * If we're here then the cache is full. The initial subword of length
             * degree - WriteCacheLetters of the current word should be the same
             * as the corresponding subword for all the tiles that we wrote to
             * the cache. Hence, it is safe to use this word as the index for
             * writing out the write cache.
             */
            cache_count = 0;
            flush_cache(word);
//            std::fill(fwd_write_cache, fwd_write_cache + cache_size, scalar_type());
//            std::fill(rev_write_cache, rev_write_cache + cache_size, scalar_type());
        }
    }

    template <unsigned D>
    void write_to_cache(const unpacked_tensor_word<Width, D>& word, index_type subtile_i, index_type subtile_j) {
        assert(word.degree() > WriteCacheLetters);
        auto ibound = base_helper::subtile_bound(subtile_i);
        auto jbound = base_helper::subtile_bound(subtile_j);

        auto degree = word.degree();
        auto index = word.split_right_index(degree - WriteCacheLetters);
        assert(index < index_type(tsi::powers[WriteCacheLetters]));

        auto offset = (subtile_i*cache_ncols + subtile_j)*tile_info::tile_width + index*tile_info::tile_stride;

        pointer ptr = fwd_write_cache + offset;
        LA_PREFETCH_ET1(ptr);
        if (tile_info::num_subtiles > 1) {
            base_helper::write_tile_assign(ptr, out_tile_ptr(), ibound, jbound, cache_ncols, tile_info::tile_width);
        } else {
            base_helper::template write_tile_assign<tile_info::tile_width, tile_info::tile_width>(ptr, out_tile_ptr(), cache_ncols);
        }

#ifdef LIBALGEBRA_TM_UPDATE_REVERSE_INLINE
        if (reverse_write_ptr != nullptr && IDEG(degree) < base::out_deg - 2*tile_letters) {
            index = word.split_right_reverse_index(degree - WriteCacheLetters);
            offset = (subtile_j*cache_ncols + subtile_i)*tile_info::tile_width + index*tile_info::tile_stride;
            ptr = rev_write_cache + offset;
            LA_PREFETCH_ET1(ptr);
            base_helper::permute_write_tile();
            if (tile_info::num_subtiles > 1) {
                base_helper::write_tile_assign(ptr, out_tile_ptr(), jbound, ibound, cache_ncols, tile_info::tile_width);
            } else {
                base_helper::template write_tile_assign<tile_info::tile_width, tile_info::tile_width>(ptr, out_tile_ptr(), cache_ncols);
            }
        }
#endif
        advance_counter(word);
    }


    LA_INLINE_ALWAYS void write_tile_impl(pointer optr, index_type stride, index_type ibound, index_type jbound)
    {
        const_pointer tptr = out_tile_ptr();

        if (base::is_inplace()) {
            base_helper::write_tile_assign(optr, tptr, ibound, jbound, stride, tile_width);
        }
        else {
            base_helper::write_tile_acc(optr, tptr, ibound, jbound, stride, tile_width);
        }
    }

    template <unsigned D>
    void write_direct(const unpacked_tensor_word<Width, D>& word, index_type subtile_i, index_type subtile_j)
    {
        auto mid_deg = word.degree();
        auto degree = mid_deg + 2*tile_info::tile_letters;

        const auto stride = tsi::powers[degree - tile_letters];
        const auto ibound = base_helper::subtile_bound(subtile_i);
        const auto jbound = base_helper::subtile_bound(subtile_j);
        const auto subtile_offset = (subtile_i * stride + subtile_j) * tile_width;

        pointer optr = base::fwd_write(degree);
        optr += word.to_index() * tile_info::tile_stride;
        optr += subtile_offset;
        LA_PREFETCH_ET1(optr);

        write_tile_impl(optr, stride, ibound, jbound);

#ifdef LIBALGEBRA_TM_UPDATE_REVERSE_INLINE
        if (reverse_write_ptr != nullptr && IDEG(degree) < base::out_deg) {
            // Write out reverse data
            optr = m_out_rev_levels[degree];
            optr += word.to_reverse_index() * tile_info::tile_stride;
            optr += (subtile_j * stride + subtile_i) * tile_width;

            LA_PREFETCH_ET1(optr);
            base_helper::permute_write_tile();

            write_tile_impl(optr, stride, jbound, ibound);
            //        base_helper::write_tile_reverse_impl(optr, stride, ibound, jbound);
            //            base_helper::write_tile(reverse_write_ptr, degree, reverse_index, subtile_i, subtile_j);
        }
#endif
    }

public:
    template <unsigned D>
    void write_tile(const unpacked_tensor_word<Width, D>& word, IDIMN subtile_i = 0, IDIMN subtile_j = 0)
    {
        if (WriteCacheLetters > 0 && word.degree() > WriteCacheLetters) {
            write_to_cache(word, subtile_i, subtile_j);
        } else {
            write_direct(word, subtile_i, subtile_j);
        }
    }

    bool cache_empty() const noexcept { return cache_count == 0; }

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


    const_pointer left_fwd_read_ptr(IDEG degree, IDIMN index, IDIMN subtile_i) const noexcept
    {
        return base::left_fwd_read(degree, index * tile_info::tile_stride + subtile_i * tile_info::tile_width /* * tsi::powers[degree-tile_letters]*/);
    }
    const_pointer left_rev_read_ptr(IDEG degree, IDIMN index, IDIMN subtile_i) const noexcept
    {
        return m_lhs_rev_levels[degree] + index * tile_info::tile_stride + subtile_i * tile_width;
    }
    const_pointer right_fwd_read_ptr(IDEG degree, IDIMN index, IDIMN subtile_j) const noexcept
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
        if (lhs.dimension() != 0 && rhs.dimension() != 0) {
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

template<DEG Width, DEG Depth, IDEG TileLetters, IDEG WriteCacheLetters>
class tiled_free_tensor_multiplication
    : public traditional_free_tensor_multiplication<Width, Depth>
{
    using base = traditional_free_tensor_multiplication<Width, Depth>;

    bool debug = false;

    template<typename C>
    using helper_type = dtl::tiled_free_tensor_multiplication_helper<
            Width, Depth, C, TileLetters == 0 ? LA_DEFAULT_TILE_PARAM(Width, Depth, typename C::S) : TileLetters, WriteCacheLetters>;


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

    template<typename C, typename Fn>
    LA_INLINE_ALWAYS static void
    linear_fwd_mul_inplace_left(pointer<C> optr, const_pointer<C> lptr, const_reference<C> runit, index_type bound, Fn op) noexcept
    {
        for (index_type i = 0; i < bound; ++i) {
            optr[i] = op(lptr[i] * runit);
        }
    }

    template<typename C, typename Fn>
    LA_INLINE_ALWAYS static void
    linear_fwd_mul_left(pointer<C> optr, const_pointer<C> lptr, const_reference<C> runit, index_type bound, Fn op) noexcept
    {
        for (index_type i = 0; i < bound; ++i) {
            optr[i] += op(lptr[i] * runit);
        }
    }

    template<typename C, typename Fn>
    LA_INLINE_ALWAYS static void
    linear_fwd_mul_inplace_right(pointer<C> optr, const_pointer<C> rptr, const_reference<C> lunit, index_type bound, Fn op) noexcept
    {
        for (index_type i = 0; i < bound; ++i) {
            optr[i] = op(lunit * rptr[i]);
        }
    }

    template<typename C, typename Fn>
    LA_INLINE_ALWAYS static void
    linear_fwd_mul_right(pointer<C> optr, const_pointer<C> rptr, const_reference<C> lunit, index_type bound, Fn op) noexcept
    {
        for (index_type i = 0; i < bound; ++i) {
            optr[i] += op(lunit * rptr[i]);
        }
    }

    template<typename C, typename Fn>
    LA_INLINE_ALWAYS static void impl_mid_compute(pointer<C> tile,
                                          const_pointer<C> lhs_tile,
                                          const_pointer<C> rhs_tile,
                                          Fn op) noexcept
    {
//        BOOST_ALIGN_ASSUME_ALIGNED(tile, LA_CACHELINE_BYTES);
        BOOST_ALIGN_ASSUME_ALIGNED(lhs_tile, LA_CACHELINE_BYTES);
        BOOST_ALIGN_ASSUME_ALIGNED(rhs_tile, LA_CACHELINE_BYTES);
        constexpr auto tile_width = helper_type<C>::tile_width;
        pointer<C> optr = tile;
        for (IDIMN i = 0; i < tile_width; ++i) {
            for (IDIMN j = 0; j < tile_width; ++j) {
                optr[j] += op(lhs_tile[i] * rhs_tile[j]);
            }
            optr += tile_width;
        }
    }

    template<typename C, typename Fn>
    LA_INLINE_ALWAYS static void impl_lb1(pointer<C> tile,
                                          const_pointer<C> lptr,
                                          const_pointer<C> rptr,
                                          index_type i1bound,
                                          index_type i2bound,
                                          Fn op)
    {
//        BOOST_ALIGN_ASSUME_ALIGNED(tile, LA_CACHELINE_BYTES);
        BOOST_ALIGN_ASSUME_ALIGNED(lptr, LA_CACHELINE_BYTES);
        BOOST_ALIGN_ASSUME_ALIGNED(rptr, LA_CACHELINE_BYTES);
        constexpr auto tile_width = helper_type<C>::tile_width;
//        auto stride = i2bound * tile_width;
        pointer<C> optr = tile;
        for (index_type i1 = 0; i1 < i1bound; ++i1) {
//            auto i = i1 * i2bound;
            for (index_type j = 0; j < tile_width; ++j) {
                optr[j] += op(lptr[i1] * rptr[j]);
            }
        }
    }

    template <typename C, index_type SBound, index_type RBound, typename Fn>
    LA_INLINE_ALWAYS static void impl_lb1(pointer<C> tptr, const_pointer<C> lptr, const_pointer<C> rptr, Fn op)
    {
        //        BOOST_ALIGN_ASSUME_ALIGNED(tile, LA_CACHELINE_BYTES);
        BOOST_ALIGN_ASSUME_ALIGNED(lptr, LA_CACHELINE_BYTES);
        BOOST_ALIGN_ASSUME_ALIGNED(rptr, LA_CACHELINE_BYTES);
        constexpr auto tile_width = helper_type<C>::tile_width;
        //        auto stride = i2bound * tile_width;
        for (index_type i1 = 0; i1 < SBound; ++i1) {
            //            auto i = i1 * i2bound;
            for (index_type j = 0; j < tile_width; ++j) {
                tptr[j] += op(lptr[i1] * rptr[j]);
            }
            tptr += RBound * tile_width;
        }
    }


    template<typename C, typename Fn>
    LA_INLINE_ALWAYS static void impl_1br(pointer<C> optr,
                                          const_pointer<C> lptr,
                                          const_pointer<C> rptr,
                                          index_type j2bound,
                                          Fn op)
    {
        constexpr auto tile_width = helper_type<C>::tile_width;
        BOOST_ALIGN_ASSUME_ALIGNED(lptr, LA_CACHELINE_BYTES);
        BOOST_ALIGN_ASSUME_ALIGNED(rptr, LA_CACHELINE_BYTES);
        for (index_type j2 = 0; j2 < j2bound; ++j2) {
            for (index_type i = 0; i < tile_width; ++i) {

                    //                auto j = j1*j2bound + j2;
                optr[i*tile_width + j2] += op(lptr[i] * rptr[j2]);
            }
//            optr += tile_width;
        }
    }

    template <typename C, index_type SBound, index_type RBound, typename Fn>
    LA_INLINE_ALWAYS static void impl_1br(pointer<C> tptr, const_pointer<C> lptr, const_pointer<C> rptr, Fn op)
    {
        constexpr auto tile_width = helper_type<C>::tile_width;
        //        BOOST_ALIGN_ASSUME_ALIGNED(tile, LA_CACHELINE_BYTES);
        BOOST_ALIGN_ASSUME_ALIGNED(lptr, LA_CACHELINE_BYTES);
        BOOST_ALIGN_ASSUME_ALIGNED(rptr, LA_CACHELINE_BYTES);
        for (index_type j2 = 0; j2 < SBound; ++j2) {
            for (index_type i = 0; i < tile_width; ++i) {
                tptr[i*tile_width + j2] += op(lptr[i] * rptr[j2]);
            }
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
        BOOST_ALIGN_ASSUME_ALIGNED(tile, LA_CACHELINE_BYTES);
        constexpr auto tile_width = helper_type<C>::tile_width;
        for (IDIMN i = 0; i < ibound; ++i) {
            for (IDIMN j = 0; j < jbound; ++j) {
                tile[i * tile_width + j] += op(lhs_fwd_ptr[i * lhs_stride] * rhs_fwd_ptr[j]);
            }
        }
    }

protected:


    template<typename Coeffs, typename Fn>
    LA_INLINE_ALWAYS static void impl_top_zipped(helper_type<Coeffs>& helper,
                                const_pointer<Coeffs> lsrc,
                                const_pointer<Coeffs> rsrc,
                                IDEG degree,
                                index_type stride,
                                index_type ibound,
                                index_type jbound,
                                Fn op)
    {
        constexpr auto tile_width = helper_type<Coeffs>::tile_width;
        pointer<Coeffs> tptr = helper.out_tile_ptr();
        const_pointer<Coeffs> lptr = helper.left_read_tile_ptr();
        const_pointer<Coeffs> rptr = helper.right_read_tile_ptr();

        const_reference<Coeffs> lunit = helper.left_unit();
        const_reference<Coeffs> runit = helper.right_unit();

//        BOOST_ALIGN_ASSUME_ALIGNED(lptr, LA_CACHELINE_BYTES);
//        BOOST_ALIGN_ASSUME_ALIGNED(rptr, LA_CACHELINE_BYTES);

        for (index_type i = 0; i < ibound; ++i) {
            LA_PREFETCH_T1(lsrc+(i+1)*stride);
            LA_PREFETCH_T1(rsrc+(i+1)*stride);
            helper.read_left_tile(lsrc + i*stride, jbound);
//            for (index_type j=0; j<jbound; ++j) {
//                tptr[i*tile_width + j] += op(lptr[j] * runit);
//            }

            helper.read_right_tile(rsrc + i*stride, jbound);
//            for (index_type j=0; j<jbound; ++j) {
//                tptr[i*tile_width + j] += op(lunit*rptr[j]);
//            }
//            lptr = lsrc + i*stride;
//            rptr = rsrc + i*stride;

            for (index_type j = 0; j < tile_width; ++j) {
                tptr[i*tile_width + j] += op(lptr[j] * runit) + op(lunit * rptr[j]);
            }
//            tptr += tile_width;
//            lsrc += stride;
//            rsrc += stride;
        }
    }

    template<typename Coeffs, typename Fn>
    static void impl_top_left_only(helper_type<Coeffs>& helper,
                                   const_pointer<Coeffs> lsrc,
                                   IDEG degree,
                                   index_type stride,
                                   index_type ibound,
                                   index_type jbound,
                                   Fn op)
    {
        constexpr auto tile_width = helper_type<Coeffs>::tile_width;
        pointer<Coeffs> tptr = helper.out_tile_ptr();
        const_pointer<Coeffs> lptr = helper.left_read_tile_ptr();


        const_reference<Coeffs> runit = helper.right_unit();

        for (index_type i = 0; i < ibound; ++i, lsrc += stride) {
            helper.read_left_tile(lsrc, jbound);

            for (index_type j = 0; j < jbound; ++j) {
                tptr[i * tile_width + j] += op(lptr[j] * runit);
            }

        }
    }

    template<typename Coeffs, typename Fn>
    static void impl_top_right_only(helper_type<Coeffs>& helper,
                                    const_pointer<Coeffs> rsrc,
                                    IDEG degree,
                                    index_type stride,
                                    index_type ibound,
                                    index_type jbound,
                                    Fn op)
    {
        constexpr auto tile_width = helper_type<Coeffs>::tile_width;
        pointer<Coeffs> tptr = helper.out_tile_ptr();
        const_pointer<Coeffs> rptr = helper.right_read_tile_ptr();

        const_reference<Coeffs> lunit = helper.left_unit();


        for (index_type i = 0; i < ibound; ++i, rsrc += stride) {
            helper.read_right_tile(rsrc, jbound);

            for (index_type j = 0; j < jbound; ++j) {
                tptr[i * tile_width + j] += op(lunit * rptr[j]);
            }
        }

    }

private:

    template <typename Coeffs, IDEG SmallCase>
    struct small_case_helper
    {
        using helper_t = helper_type<Coeffs>;
        using next_t = small_case_helper<Coeffs, SmallCase-1>;
        static_assert(0<SmallCase && SmallCase <= helper_t::tile_letters,
                      "Small cases must be of depth < tile_letters");

        static constexpr index_type small_bound = integer_maths::power(index_type(Width), SmallCase);
        static constexpr index_type remainder_bound = integer_maths::power(index_type(Width), helper_t::tile_letters - SmallCase);


        template <typename Fn>
        static void left(helper_t& helper, IDEG out_deg, IDEG lhs_deg, index_type k, Fn op)
        {
            if (lhs_deg == SmallCase) {
                const auto rhs_deg = out_deg - lhs_deg;
                pointer<Coeffs> tptr = helper.out_tile_ptr();
                const_pointer<Coeffs> lptr = helper.left_read_tile_ptr();
                const_pointer<Coeffs> rptr = helper.right_read_tile_ptr();

                helper.template read_left_tile<remainder_bound>(helper.left_fwd_read(lhs_deg));
                const_pointer<Coeffs> rhs_read = helper.right_fwd_read_ptr(rhs_deg, k, 0);

                const auto key_offset = static_cast<index_type>(tsi::powers[out_deg - 2 * helper_t::tile_letters]);

//                const_pointer<Coeffs>* reads = helper.read_queue();
//                for (index_type i2 = 0; i2 < remainder_bound; ++i2) {
//                    reads[i2] = rhs_read + i2 * key_offset * helper_t::tile_stride;
//                    LA_PREFETCH_T1(reads[i2]);
//                }

                for (index_type i2 = 0; i2 < remainder_bound; ++i2) {
                    //            helper.read_right_tile(rhs_read + i2*key_offset);
                    helper.template read_right_tile<helper_t::tile_width>(rhs_read + i2*key_offset*helper_t::tile_stride);
//                    rhs_read = helper.right_fwd_read_ptr(rhs_deg, (i2+1)*key_offset + k, 0);
                    LA_PREFETCH_T1(rhs_read + (i2+1)*key_offset);
//                    helper.read_right_tile(reads[i2]);
//                    impl_lb1<Coeffs>(tptr + i2*helper_type<Coeffs>::tile_width, lptr, rptr, small_bound, remainder_bound, op);
                    impl_lb1<Coeffs, small_bound, remainder_bound>(tptr + i2 * helper_t::tile_width, lptr, rptr, op);
                }
            }
            else {
                next_t::left(helper, out_deg, lhs_deg, k, op);
            }
        }

        template<typename Fn>
        static void right(helper_t& helper, IDEG out_deg, IDEG rhs_deg, index_type k_reverse, Fn op)
        {
            if (rhs_deg == SmallCase) {
                const auto lhs_deg = out_deg - rhs_deg;
                const auto mid_deg = out_deg - 2 * helper_t::tile_letters;

                helper.template read_right_tile<remainder_bound>(helper.right_fwd_read(rhs_deg));

                pointer<Coeffs> tptr = helper.out_tile_ptr();
                const_pointer<Coeffs> lptr = helper.left_read_tile_ptr();
                const_pointer<Coeffs> rptr = helper.right_read_tile_ptr();

                const auto key_offset = static_cast<index_type>(tsi::powers[mid_deg]);

//                const_pointer<Coeffs>* reads = helper.read_queue();
//                const_pointer<Coeffs> reads[remainder_bound];
//                const_pointer<Coeffs> left_read = helper.left_rev_read_ptr(lhs_deg, k_reverse, 0);
                unpacked_tensor_word<Width, helper_t::tile_letters> j1_word;
                j1_word.reset(helper_t::tile_letters - SmallCase);
//                for (index_type j1 = 0; j1 < remainder_bound; ++j1, ++j1_word) {
//                    auto rj1 = j1_word.to_reverse_index();
//                    reads[j1] = left_read + rj1 * key_offset * helper_t::tile_stride;
//                    LA_PREFETCH_T1(reads[j1]);
//                }

                index_type rj1 = 0;
                const_pointer<Coeffs> left_read = helper.left_rev_read_ptr(lhs_deg, k_reverse, 0);
                LA_PREFETCH_T1(left_read);
                for (index_type j1 = 0; j1 < remainder_bound; ++j1) {
//                                auto rj1 = helper.reverse_key(j1deg, j1);
                    helper.template read_left_tile<helper_t::tile_width>(left_read);
                    ++j1_word;
                    rj1 = j1_word.to_reverse_index();
                    left_read = helper.left_rev_read_ptr(lhs_deg, rj1*key_offset + k_reverse, 0);
                    LA_PREFETCH_T1(left_read);
                    //                                helper.read_left_tile(left_read + rj1*key_offset);
//                    helper.read_left_tile(reads[j1]);
                    helper.permute_left_tile();
//                    impl_1br<Coeffs>(tptr + j1*small_bound, lptr, rptr, small_bound, op);
                    impl_1br<Coeffs, small_bound, remainder_bound>(tptr + j1 * small_bound, lptr, rptr, op);
                }
            } else {
                next_t::right(helper, out_deg, rhs_deg, k_reverse, op);
            }
        }

    };

   template <typename Coeffs>
   struct small_case_helper<Coeffs, 0>
   {
        using helper_t = helper_type<Coeffs>;

        template<typename Fn>
        static void left(helper_t& helper, IDEG out_deg, IDEG lhs_deg, index_type k, Fn op)
        {}

        template<typename Fn>
        static void right(helper_t& helper, IDEG out_deg, IDEG rhs_deg, index_type k_reverse, Fn op)
        {}
   };


public:

    template<typename Coeffs, typename Fn>
    static void
    impl_lhs_small(helper_type<Coeffs>& helper,
                   IDEG out_deg,
                   IDEG lhs_deg,
                   index_type k,
                   Fn op) noexcept
    {
        small_case_helper<Coeffs, helper_type<Coeffs>::tile_letters>::left(helper, out_deg, lhs_deg, k, op);


//        const auto rhs_deg = out_deg - lhs_deg;
//        constexpr auto tile_letters = helper_type<Coeffs>::tile_letters;
//
//        assert(1 <= lhs_deg && lhs_deg <= tile_letters);
//
//        const auto i1bound = static_cast<index_type>(tsi::powers[tile_letters - lhs_deg]);
//        const auto i2bound = static_cast<index_type>(tsi::powers[lhs_deg]);
//
//        pointer<Coeffs> tptr = helper.out_tile_ptr();
//        const_pointer<Coeffs> lptr = helper.left_read_tile_ptr();
//        const_pointer<Coeffs> rptr = helper.right_read_tile_ptr();
//        helper.read_left_tile(helper.left_fwd_read(lhs_deg), i1bound);
//        const_pointer<Coeffs> rhs_read = helper.right_fwd_read_ptr(rhs_deg, k, 0);
//
//        const auto key_offset = static_cast<index_type>(tsi::powers[out_deg - 2*tile_letters]);
//
//        const_pointer<Coeffs>* reads = helper.read_queue();
//        for (index_type i2=0; i2<i2bound; ++i2) {
//            reads[i2] = helper.right_fwd_read_ptr(rhs_deg, i2*key_offset + k, 0);
//            LA_PREFETCH_T1(reads[i2]);
//        }
//
//
//        for (index_type i2 = 0; i2 < i2bound; ++i2) {
////            helper.read_right_tile(rhs_read + i2*key_offset);
//            helper.read_right_tile(reads[i2]);
//            impl_lb1<Coeffs>(tptr + i2*helper_type<Coeffs>::tile_width, lptr, rptr, i1bound, i2bound, op);
//        }
    }

    template <typename Coeffs, typename Array, typename Fn>
    static void impl_mid(helper_type<Coeffs>& helper,
                         IDEG out_deg,
                         IDEG lhs_deg,
                         const Array& left_reads,
                         const Array& right_reads,
                         index_type ibound, index_type jbound,
                         Fn fn)
    {
        helper.read_left_tile(left_reads[lhs_deg], ibound);
        helper.read_right_tile(right_reads[out_deg - lhs_deg], jbound);
        helper.permute_left_tile();
        impl_mid_compute<Coeffs>(helper.out_tile_ptr(), helper.left_read_tile_ptr(), helper.right_read_tile_ptr(), fn);
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
    impl_rhs_small(helper_type<Coeffs>& helper,
                           IDEG out_deg,
                           IDEG lhs_deg,
                           IDIMN k_reverse,
                           Fn op) noexcept
    {
        small_case_helper<Coeffs, helper_type<Coeffs>::tile_letters>::right(helper, out_deg, out_deg - lhs_deg, k_reverse, op);
//        constexpr auto tile_letters = helper_type<Coeffs>::tile_letters;
//
//        const auto rhs_deg = out_deg - lhs_deg;
//
//        const auto mid_deg = out_deg - 2 * tile_letters;
//        const auto j1deg = tile_letters - rhs_deg;
//        const auto j1bound = static_cast<index_type>(tsi::powers[j1deg]);
//        const auto j2bound = static_cast<index_type>(tsi::powers[rhs_deg]);
//        helper.read_right_tile(helper.right_fwd_read(rhs_deg), j2bound);
//
//        pointer<Coeffs> tptr = helper.out_tile_ptr();
//        const_pointer<Coeffs> lptr = helper.left_read_tile_ptr();
//        const_pointer<Coeffs> rptr = helper.right_read_tile_ptr();
//
//        const auto key_offset = static_cast<index_type>(tsi::powers[mid_deg]);
//
//        const_pointer<Coeffs>* reads = helper.read_queue();
//        const_pointer<Coeffs> left_read = helper.left_rev_read_ptr(lhs_deg, k_reverse, 0);
//        unpacked_tensor_word<Width, tile_letters> j1_word;
//        j1_word.reset(j1deg);
////        for (index_type j1=0; j1<j1bound; ++j1, ++j1_word) {
////            auto rj1 = j1_word.to_reverse_index();
////            reads[j1] = helper.left_rev_read_ptr(lhs_deg, rj1 * key_offset + k_reverse, 0);
////            LA_PREFETCH_T1(reads[j1]);
////        }
//
//        index_type rj1 = 0;
//        for (index_type j1 = 0; j1 < j1bound; ++j1) {
////            auto rj1 = helper.reverse_key(j1deg, j1);
////            helper.template read_left_tile<helper_type<Coeffs>::tile_width>(reads[j1]);
//            const_pointer<Coeffs> left_read = helper.left_rev_read_ptr(lhs_deg, rj1*key_offset + k_reverse, 0);
//            helper.read_left_tile(left_read);
//            ++j1_word;
//            rj1 = j1_word.to_reverse_index();
//            helper.permute_left_tile();
//            impl_1br<Coeffs>(tptr, lptr, rptr, j2bound, op);
//            tptr += j2bound;
//        }

    }

    template<typename Coeffs, typename Fn>
    static void
    impl_rhs_small_no_reverse(helper_type<Coeffs>& helper,
                              IDEG out_deg,
                              IDEG lhs_deg,
                              IDIMN k,
                              Fn op) noexcept
    {
        constexpr auto tile_width = helper_type<Coeffs>::tile_width;
        constexpr auto tile_letters = helper_type<Coeffs>::tile_letters;

        const auto rhs_deg = out_deg - lhs_deg;
        const auto split_left_letters = tile_letters - rhs_deg;
        const auto lhs_stride = tsi::powers[lhs_deg - tile_letters];

        for (IDIMN j = 0; j < tile_width; ++j) {
            const auto split = helper.split_key(rhs_deg, j);

            const auto& right_val = *helper.right_fwd_read(rhs_deg, split.second);
            const auto lhs_key = helper.combine_keys(split_left_letters, k, split.first);

            impl_ulmd<Coeffs>(helper.out_tile_ptr(),
                              helper.left_fwd_read(lhs_deg, lhs_key),
                              right_val,
                              op,
                              j,
                              tile_width,
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

        //        impl_outer_cases(helper, op);

        unpacked_tensor_word<Width, Depth> word;
        for (IDEG out_deg = max_degree; out_deg > 2 * tile_letters; --out_deg) {
            const auto mid_deg = out_deg - 2 * tile_letters;
            const auto mid_end = out_deg - tile_letters;
            const auto stride = static_cast<IDIMN>(tsi::powers[out_deg - tile_letters]);

            auto lhs_deg_min = std::max(IDEG(1), out_deg - rhs_max_deg);
            auto lhs_deg_max = std::min(out_deg - 1, old_lhs_deg);

            word.reset(mid_deg);

            for (IDIMN k = 0; k < static_cast<IDIMN>(tsi::powers[mid_deg]); ++k, ++word) {
                assert(k == word.to_index());
                assert(IDEG(word.degree()) == mid_deg);
                auto k_reverse = word.to_reverse_index();

                for (IDIMN subtile_i = 0; subtile_i < num_subtiles; ++subtile_i) {
                    const auto ibound = helper.subtile_bound(subtile_i);
                    for (IDIMN subtile_j = 0; subtile_j < num_subtiles; ++subtile_j) {
                        const auto jbound = helper.subtile_bound(subtile_j);
                        const auto subtile_offset = (subtile_i*stride + subtile_j)*tile_width;

//                        _mm_prefetch(helper.out_tile_ptr(), _MM_HINT_T0);
//                        _mm_prefetch(helper.left_read_tile_ptr(), _MM_HINT_T0);
//                        _mm_prefetch(helper.right_read_tile_ptr(), _MM_HINT_T0);

                        helper.reset_tile(out_deg, k, k_reverse, subtile_i, subtile_j);

                        // Setup read pointers
                        const auto& rhs_unit = *right_reads[0];
                        const auto& lhs_unit = *left_reads[0];
                        bool left_ok = out_deg <= old_lhs_deg && rhs_unit != Coeffs::zero;
                        bool right_ok = out_deg <= rhs_max_deg && lhs_unit != Coeffs::zero;                       //                        if (out_deg < max_degree) {
                        if (left_ok) {
                            left_reads[out_deg] = helper.left_fwd_read(out_deg, k * helper.tile_stride + subtile_offset);
                        }
                        if (right_ok) {
                            right_reads[out_deg] = helper.right_fwd_read(out_deg, k * helper.tile_stride + subtile_offset);
                        }
                        //                        }
                        for (IDEG i = std::max(tile_letters, lhs_deg_min); i <= std::min(out_deg - tile_letters, lhs_deg_max); ++i) {
//                            auto split = helper.split_key(out_deg - i - tile_letters, k);
//                            auto lkey = helper.reverse_key(i - tile_letters, split.first);
//                            auto rkey = split.second;
                            auto lkey = word.split_left_reverse_index(i - tile_letters);
                            auto rkey = word.split_right_index(i-tile_letters);
                            left_reads[i] = helper.left_rev_read_ptr(i, lkey, subtile_i);
                            right_reads[out_deg - i] = helper.right_fwd_read_ptr(out_deg - i, rkey, subtile_j);
                        }

                        //                        if (out_deg < max_degree) {



                        for (IDEG lhs_deg = lhs_deg_min; lhs_deg <= lhs_deg_max; ++lhs_deg) {
                            if (lhs_deg < tile_letters) {
                                impl_lhs_small(helper, out_deg, lhs_deg, k, op);
                            }
                            else if (lhs_deg <= mid_end && lhs_deg < old_lhs_deg) {
                                impl_mid(helper, out_deg, lhs_deg, left_reads, right_reads, ibound, jbound, op);
                                //                                impl_mid_cases_reverse(helper, op, out_deg, lhs_deg, k, subtile_i, subtile_j);
                            }
                            else if (lhs_deg <= mid_end && lhs_deg == old_lhs_deg) {
                                impl_mid_cases_no_reverse(helper, op, out_deg, lhs_deg, k, subtile_i, subtile_j);
                            }
                            else if (lhs_deg > mid_end && lhs_deg < old_lhs_deg) {
                                impl_rhs_small(helper, out_deg, lhs_deg, k_reverse, op);
                            }
                            else if (lhs_deg > mid_end) {
                                impl_rhs_small_no_reverse(helper, out_deg, lhs_deg, k, op);
                            }
                        }

                        if (left_ok && right_ok) {
                            impl_top_zipped(helper, left_reads[out_deg], right_reads[out_deg],
                                            out_deg, stride, ibound, jbound, op);
                        }
                        else if (left_ok) {
                            impl_top_left_only(helper, left_reads[out_deg],
                                               out_deg, stride, ibound, jbound, op);
                        }
                        else if (right_ok) {
                            impl_top_right_only(helper, right_reads[out_deg],
                                                out_deg, stride, ibound, jbound, op);
                        }

                        helper.write_tile(word, subtile_i, subtile_j);
                    }// subtile_j
                }    // subtile_i
            }        // k
            assert(helper.cache_empty());

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

    template<typename B1, typename Coeffs, typename Fn, typename OriginalVectors>
    void fma(vectors::dense_vector<B1, Coeffs>& out,
             const vectors::dense_vector<B1, Coeffs>& lhs,
             const vectors::dense_vector<B1, Coeffs>& rhs,
             Fn op,
             DEG max_degree,
             OriginalVectors& orig) const
    {
        if (max_degree <= 2 * helper_type<Coeffs>::tile_letters || lhs.degree() + rhs.degree() <= 2 * helper_type<Coeffs>::tile_letters) {
            base::fma(out, lhs, rhs, op, max_degree, orig);
            return;
        }

        if (lhs.dimension() != 0 && rhs.dimension() != 0) {
            helper_type<Coeffs> helper(out, lhs, rhs, max_degree);
            fma_impl(helper, op, helper.out_degree());
            if (helper.out_degree() > 0) {
#ifdef LIBALGEBRA_TM_UPDATE_REVERSE_INLINE
                auto update_deg = std::min(helper.out_degree() - 1, 2*helper.tile_letters);
#else
                auto update_deg = helper.out_degree() - 1;
#endif
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
        if (max_degree <= 2 * helper_type<Coeffs>::tile_letters || lhs.degree() + rhs.degree() <= 2 * helper_type<Coeffs>::tile_letters) {
            base::multiply_inplace(lhs, rhs, op, max_degree, orig);
            return;
        }

        if (rhs.dimension() != 0) {
            helper_type<Coeffs> helper(lhs, rhs, max_degree);
            multiply_inplace_impl(helper, op, helper.out_degree());
            if (helper.out_degree() > 0) {
#ifdef LIBALGEBRA_TM_UPDATE_REVERSE_INLINE
                auto update_deg = std::min(helper.out_degree() - 1, 2*helper.tile_letters);
#else
                auto update_deg = helper.out_degree() - 1;
#endif
                base::update_reverse_data(lhs, update_deg);
            }
        }
        else {
            lhs.clear();
        }
    }
};


template <DEG Depth, IDEG TileLetters, IDEG WriteCacheLetters>
class tiled_free_tensor_multiplication<1, Depth, TileLetters, WriteCacheLetters>
        : public traditional_free_tensor_multiplication<1, Depth>
{};


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
    result_type operator()(argument_type lhs, argument_type rhs) const
    {
        static const boost::container::small_vector<pair_type, 0> null;

        if ((lhs.size() + rhs.size()) > Depth) {
            return null;
        }

        return base::cached_compute(lhs, rhs);
        //        return half_shuffle_base::shuffle(lhs, rhs);
    }
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
    result_type operator()(argument_type lhs, argument_type rhs) const
    {
        static const boost::container::small_vector<pair_type, 0> null;

        if (Depth > 0 && (lhs.size() + rhs.size()) > Depth) {
            return null;
        }

        return base::cached_compute(lhs, rhs);
    }
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

template<typename Coeff, DEG n_letters, DEG max_degree,
         template <typename, typename, typename...> class VectorType,
        typename...>
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
         template <DEG, DEG> class FTMultiplication,
         typename... Args>
class free_tensor
    : public algebra<
              free_tensor_basis<n_letters, max_degree>,
              Coeff,
              FTMultiplication<n_letters, max_degree>,
              VectorType,
              free_tensor<Coeff, n_letters, max_degree, VectorType, FTMultiplication, Args...>,
              Args...>
{
    typedef FTMultiplication<n_letters, max_degree> multiplication_t;

    using base = algebra<
            free_tensor_basis<n_letters, max_degree>,
            Coeff,
            FTMultiplication<n_letters, max_degree>,
            VectorType,
            free_tensor,
            Args...>;

    template<template<typename, typename, typename...> class VT, template <DEG, DEG> class VFTM, typename... VArgs>
    static void resize_for_degree(free_tensor<Coeff, n_letters, max_degree, VT, VFTM, VArgs...>& arg, DEG degree)
    {}

    template<typename... VArgs>
    static void resize_for_degree(free_tensor<Coeff, n_letters, max_degree, ::alg::vectors::dense_vector, FTMultiplication, VArgs...>& arg, DEG degree)
    {
        arg.base_vector().resize_to_degree(degree);
    }

public:
    /// The basis type.
    typedef free_tensor_basis<n_letters, max_degree> BASIS;
    /// Import of the KEY type.
    typedef typename BASIS::KEY KEY;
    /// The algebra type.
    typedef algebra<BASIS, Coeff, multiplication_t, VectorType, free_tensor<Coeff, n_letters, max_degree, VectorType, FTMultiplication, Args...>, Args...> ALG;

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

    free_tensor() : base() {}

    explicit free_tensor(typename boost::call_traits<scalar_type>::param_type s)
        : base(key_type{}, s)
    {}

    template<typename Letter, typename Scalar>
    explicit free_tensor(Letter let, Scalar sca)
        : base(base::basis.keyofletter(LET(let)), scalar_type(sca))
    {}

    /// Computes the truncated exponential of a free_tensor instance.

    template <typename Tensor>
    inline friend typename std::enable_if<std::is_base_of<free_tensor, Tensor>::value, Tensor>::type
     exp(const Tensor& arg)
    {
        // Computes the truncated exponential of arg
        // 1 + arg + arg^2/2! + ... + arg^n/n! where n = max_degree
        KEY kunit;
        Tensor tunit(kunit);
        Tensor result(tunit);

        resize_for_degree(result, max_degree);

        for (DEG i = max_degree; i >= 1; --i) {
            result.mul_scal_div(arg, typename Coeff::Q(i));
            result += tunit;
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
    template <typename Tensor>
    inline friend typename std::enable_if<std::is_base_of<free_tensor, Tensor>::value, Tensor>::type
    log(const Tensor& arg)
    {
        // Computes the truncated log of arg up to degree max_degree
        // The coef. of the constant term (empty word in the monoid) of arg
        // is forced to 1.
        // log(arg) = log(1+x) = x - x^2/2 + ... + (-1)^(n+1) x^n/n.
        // max_degree must be > 0

        KEY kunit;
        Tensor tunit(kunit);
        Tensor x(arg);
        iterator it = x.find(kunit);
        if (it != x.end()) {
            x.erase(it);
        }
        Tensor result;

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
    template<typename Tensor>
    friend typename std::enable_if<std::is_base_of<free_tensor, Tensor>::value, Tensor>::type
    inverse(const Tensor& arg)
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
        SCA a(0);
        Tensor x, z(a);

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
        Tensor free_tensor_a_inverse(SCA(1) / a), result(free_tensor_a_inverse);
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
template<typename Coeff, DEG n_letters, DEG max_degree,
         template <typename, typename, typename...> class VectorType,
         typename... ExtraArgs>
class shuffle_tensor : public algebra<
                               shuffle_tensor_basis<n_letters, max_degree>,
                               Coeff,
                               shuffle_tensor_multiplication<n_letters, max_degree>,
                               VectorType,
                               shuffle_tensor<Coeff, n_letters, max_degree, VectorType, ExtraArgs...>
                               >
{
    typedef shuffle_tensor_multiplication<n_letters, max_degree> multiplication_t;

public:
    /// The basis type.
    typedef shuffle_tensor_basis<n_letters, max_degree> BASIS;
    /// Import of the KEY type.
    typedef typename BASIS::KEY KEY;
    /// The algebra type.
    typedef algebra<BASIS, Coeff, multiplication_t, VectorType, shuffle_tensor> ALG;

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
    template<template<typename, typename, typename...> class OVectorType, template <DEG, DEG> class FTM, typename... Args>
    shuffle_tensor(const free_tensor<Coeff, n_letters, max_degree, OVectorType, FTM, Args...>& t)
    {
        typename free_tensor<Coeff, n_letters, max_degree, OVectorType, FTM, Args...>::const_iterator i;
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

    key_type key() const
    {
        assert(index() < p_vector->dimension());
        return basis_type::index_to_key(index());
    }

    reference value() const
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

    key_type key() const
    {
        assert(index() < p_vector->dimension());
        return basis_type::index_to_key(index());
    }

    reference value() const
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
