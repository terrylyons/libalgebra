#include <vector>

#include <libalgebra/vectors/vectors.h>

namespace alg {

    namespace dtl {

        template<typename Coeff, DEG Width, DEG MaxDepth, DEG TileLetters>
        class GilesMultiplier {
            using S = typename Coeff::SCA;

            std::vector<S> left_read_tile, right_read_tile;
            std::vector<S> output_tile;
            std::vector<S> lhs_reverse_data;
            DEG m_lhs_deg, m_rhs_deg;

            using tensor_basis_type = free_tensor_basis<Width, MaxDepth>;

            using tensor_type = free_tensor <Coeff, Width, MaxDepth, vectors::dense_vector>;

        public:
            static const std::vector<DIMN> powers;

            S *output_ptr;
            const S *left_forward_read_ptr;
            const S *left_reverse_read_ptr;
            const S *right_forward_read_ptr;

            static constexpr DIMN tile_width = integer_maths::power(Width, TileLetters);
            static constexpr DIMN tile_size = tile_width * tile_width;
            static constexpr DIMN tile_shift = integer_maths::power(Width, TileLetters - 1);

            GilesMultiplier(tensor_type &out_tensor, const tensor_type &lhs_tensor, const tensor_type &rhs_tensor)
                    : left_read_tile(tile_width),
                      m_lhs_deg(lhs_tensor.degree()),
                      m_rhs_deg(rhs_tensor.degree()),
                      right_read_tile(tile_width),
                      output_tile(tile_size),
                      lhs_reverse_data(tensor_basis_type::start_of_degree(MaxDepth))
            {
                vectors::dtl::vector_base_access::convert(out_tensor).resize_to_dimension(tensor_basis_type::start_of_degree(MaxDepth+1));
                using traits = vectors::dtl::data_access<vectors::dense_vector<tensor_basis_type, Coeff>>;
                output_ptr = traits::range_begin(vectors::dtl::vector_base_access::convert(out_tensor));
                left_forward_read_ptr = traits::range_begin(vectors::dtl::vector_base_access::convert(lhs_tensor));
                right_forward_read_ptr = traits::range_begin(vectors::dtl::vector_base_access::convert(rhs_tensor));
                left_reverse_read_ptr = lhs_reverse_data.data();


                /*
                 * This logic constructs a temporary reverse data within the struct so we can access it in the algorithm
                 * Later we can replace this with logic to use the lhs_tensor's reverse data if it has some.
                 */
#ifdef LIBALGEBRA_MAX_TILE_LETTERS
                constexpr DEG CalcLetters = integer_maths::logN(static_cast<unsigned>(LIBALGEBRA_L1_CACHE_SIZE) / sizeof(S), Width) / 2;
            constexpr DEG BlockLetters = (CalcLetters > LIBALGEBRA_MAX_TILE_LETTERS) ? LIBALGEBRA_MAX_TILE_LETTERS : CalcLetters;
#else
                constexpr DEG BlockLetters = integer_maths::logN(LIBALGEBRA_L1_CACHE_SIZE / sizeof(S), Width) / 2;
#endif

                dtl::tiled_inverse_operator<Width, MaxDepth - 1, BlockLetters, S, dtl::non_signing_signer> reverser;
                reverser(left_forward_read_ptr, lhs_reverse_data.data(), lhs_tensor.degree()-1);

            }

            DEG lhs_degree() const noexcept { return m_lhs_deg; }
            DEG rhs_degree() const noexcept { return m_rhs_deg; }

            S* out_tile_ptr() noexcept
            {
                return output_tile.data();
            }
            const S* left_read_tile_ptr() const noexcept
            {
                return left_read_tile.data();
            }
            const S* right_read_tile_ptr() const noexcept
            {
                return right_read_tile.data();
            }

            void read_left_tile(DEG degree, DIMN index) noexcept
            {
                const auto start_of_degree = tensor_basis_type::start_of_degree(degree);
                const auto* ptr_begin = left_reverse_read_ptr + index * tile_width + start_of_degree;
                std::copy(ptr_begin, ptr_begin + tile_width, left_read_tile.data());
            }
            void read_right_tile(DEG degree, DIMN index) noexcept
            {
                const auto start_of_degree = tensor_basis_type::start_of_degree(degree);
                const auto* ptr_begin = right_forward_read_ptr + index * tile_width + start_of_degree;
                std::copy(ptr_begin, ptr_begin + tile_width, right_read_tile.data());
            }
            void write_tile(DEG degree, DIMN index, DIMN reverse_index) noexcept
            {
                const auto start_of_degree = tensor_basis_type::start_of_degree(degree);
                auto* optr = output_ptr + index * tile_width + start_of_degree;
                const auto* tptr = output_tile.data();
                auto stride = powers[degree-TileLetters];

                assert((tile_width-1)*stride + index*tile_width  + tile_width-1 < powers[degree]);
                for (DIMN i = 0; i < tile_width; ++i) {
                    for (DIMN j = 0; j < tile_width; ++j) {
                        optr[i*stride+j] += tptr[i*tile_width+j];
                    }
                }

                if (degree < MaxDepth) {
                    // Write out to the reverse data too
                }
            }

            static std::pair<DIMN, DIMN> split_key(DEG split_degree, DIMN key) noexcept
            {
                const DIMN splitter = powers[split_degree];
                return {key / splitter, key % splitter};
            }
            static DIMN combine_keys(DEG right_degree, DIMN left, DIMN right) noexcept
            {
                const DIMN shift = powers[right_degree];
                return left * shift + right;
            }

            static DIMN reverse_key(DEG degree, DIMN index) noexcept
            {
                if (degree < 2) {
                    return index;
                }
                else if (degree == 2) {
                    return (index % Width) * Width + (index / Width);
                }
                else {
                    assert(degree < powers.size());
                    auto left_letter = index / powers[degree - 1];
                    auto right_letter = index % Width;
                    auto middle_letters = (index / Width);
                    return combine_keys(1, combine_keys(degree - 1, right_letter, reverse_key(degree - 2, middle_letters)), left_letter);
                }
            }

            static constexpr DIMN reverse(DIMN index) noexcept
            {
                return dtl::reversing_permutation<Width, TileLetters>::permute_idx(index);
            }

            const S& left_unit() const noexcept { return left_forward_read_ptr[0]; }
            const S& right_unit() const noexcept { return right_forward_read_ptr[0]; }
            const S* left_fwd_read(DEG degree, DIMN index) const noexcept
            {
                auto offset = tensor_basis_type::start_of_degree(degree);
                return left_forward_read_ptr + index*tile_width + offset;
            }
            const S* right_fwd_read(DEG degree, DIMN index) const noexcept
            {
                const DIMN offset = tensor_basis_type::start_of_degree(degree);
                return right_forward_read_ptr + index*tile_width + offset;
            }
            const S* left_reverse_read(DEG degree, DIMN index) const noexcept
            {
                const DIMN offset = tensor_basis_type::start_of_degree(degree);
                return left_reverse_read_ptr + index + offset;
            }
            S* fwd_write(int degree) const noexcept
            {
                return output_ptr + tensor_basis_type::start_of_degree(static_cast<DEG>(degree));
            }

        }; // GilesMultiplier class

        template<typename S, DEG Width, DEG MaxDepth, DEG TileLetters>
        const std::vector<DIMN> GilesMultiplier<S, Width, MaxDepth, TileLetters>::powers = []() {
            std::vector<DIMN> result;
            result.reserve(MaxDepth + 1);
            result.push_back(1);

            for (DEG i = 1; i <= MaxDepth; ++i) {
                result.push_back(Width * result.back());
            }

            return result;
        }();

    } //namespace dtl
} // namespace alg
