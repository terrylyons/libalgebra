/* *************************************************************

Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurk� and Arend Janssen.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */

//  tensor.h

// Include once wrapper
#ifndef DJC_COROPA_LIBALGEBRA_TENSORH_SEEN
#define DJC_COROPA_LIBALGEBRA_TENSORH_SEEN

#ifndef LIBALGEBRA_L1_CACHE_SIZE
#define LIBALGEBRA_L1_CACHE_SIZE 32768   // 32Kb should be fairly standard
#endif

//#include <omp.h>

#include <unordered_set>
#include <algorithm>

#include <libalgebra/vectors/base_vector.h>

namespace alg {

namespace dtl{

    template <unsigned Width, unsigned Level>
    struct reversing_permutation
    {
        using size_type  = size_t;
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
                for (auto& v: w) {
                    result *= Width;
                    result += v;
                }
                return result;
            };

            result.reserve(tile_size); // over-sized, but that's fine
            std::unordered_set<size_type> seen;

            for (size_type lhs = 0; lhs<tile_size; ++lhs) {

                size_type rhs = idx(rword);
                if (lhs!=rhs && seen.count(lhs)==0 && seen.count(rhs)==0) {
                    result.emplace_back(lhs, rhs);
                }

                seen.insert(lhs);
                seen.insert(rhs);

                for (int i = 0; i<Level; ++i) {
                    if (word[Level-i-1]<Width-1) {
                        ++word[Level-i-1];
                        ++rword[i];
                        break;
                    }
                    else {
                        word[Level-i-1] = 0;
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
            static_assert(Level-2 > 0, "Level must be at least 3 in this specialisation");
            using next = reversing_permutation<Width, Level-2>;

            constexpr size_type shift = power(Width, Level-1);
            return last_letter(idx)*shift + next::permute_idx(middle_word(idx)) * Width + first_letter(idx);
        }


        template<typename T>
        void operator()(T* __restrict tile)  const noexcept
        {
            T tmp;

            for (size_type i=0; i < power(Width, Level); ++i)
            {
                auto j = permute_idx(i);
                if (j > i)
                {
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
        void operator()(T * __restrict tile) const noexcept
        {
            T tmp;

            for (size_type i = 0; i < power(Width, Level); ++i) {
                auto j = permute_idx(i);
                if (j > i)
                {
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

    template <DEG Width, DEG MaxDepth, DEG BlockLetters, typename Scalar>
    class tiled_inverse_operator{

        typedef Scalar SCA;

        using BASIS = free_tensor_basis<Width, MaxDepth>;

    public:
        static constexpr DEG block_letters = BlockLetters;
        static constexpr size_t block_width = power(Width, BlockLetters);
        static constexpr size_t block_size = power(Width, 2*BlockLetters);
        static constexpr unsigned max_middle_word_length = MaxDepth - 2*BlockLetters;
//        static constexpr size_type middle_word_count = tensor_alg_size(max_middle_word_length);
        static constexpr size_t block_offset = power(Width, BlockLetters);

        template <DEG Level, DEG MaxLevel>
        struct recursive_untiled_compute
        {
            using permutation_t = reversing_permutation<Width, Level>;
//            using signer_t = signing_operator<Level % 2>;
            using next_t = recursive_untiled_compute<Level+1, MaxLevel>;
            static constexpr size_t level_size = power(Width, Level);

            /*
             * Everything here is supposed to fit in cache, so we really don't need to worry
             * about locality etc when manipulating arrays. This should be very fast.
             */
            void operator()(const SCA* __restrict src_ptr, SCA* __restrict dst_ptr, DEG curr_deg) const noexcept
            {
                // Copy from src to test and adjust sign.
//                signer_t signer;
                for (size_t i=0; i<level_size; ++i) {
                    if (Level % 2 == 0)
                    {
                        dst_ptr[i] = src_ptr[i];
                        // dst_ptr[i] = signer(src_ptr[i]);
                    }
                    else
                    {
                        dst_ptr[i] = -src_ptr[i];
                        // dst_ptr[i] = signer(src_ptr[i]);
                    }
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

        template <DEG MaxLevel>
        struct recursive_untiled_compute<MaxLevel, MaxLevel>
        {
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
                for (size_t i = 0; i<level_size; ++i) {
                    if (Level % 2 == 0)
                    {
                        dst_ptr[i] = src_ptr[i];
                        // dst_ptr[i] = signer(src_ptr[i]);
                    }
                    else
                    {
                        dst_ptr[i] = -src_ptr[i];
                        // dst_ptr[i] = signer(src_ptr[i]);
                    }
                }

                // Operate on the pointer as if it were a tile of size Width^Level
                permutation_t permutation;
                permutation(dst_ptr);


                // No more recursing!
            }
        };

        recursive_untiled_compute<0U, 2*BlockLetters-1> recurse;

        void operator()(const SCA* src_ptr, SCA* dst_ptr, const unsigned curr_degree) const noexcept
        {

            if (src_ptr == nullptr) // if pointer to source is null
            {
                return;
            }

            recurse(src_ptr, dst_ptr, curr_degree);

            for (unsigned int length = 0; length <= max_middle_word_length && length + 2*block_letters <= curr_degree ; ++length)
            {
                auto istart = BASIS::start_of_degree(length);
                auto iend = BASIS::start_of_degree(length + 1);

                auto src_dst_offset = BASIS::start_of_degree(length + 2 * BlockLetters);

                auto src_p = src_ptr + src_dst_offset;
                auto dst_p = dst_ptr + src_dst_offset;

                auto key_start = BASIS::index_to_key(istart);
                auto key_end = BASIS::index_to_key(iend);

                auto word_idx = istart;

//#pragma omp parallel for
                for (auto word = key_start; word != key_end; word = BASIS::nextkey(word), ++word_idx)
                {
                    auto rword_index = BASIS::key_to_index(word.reverse());

                    if (length % 2 == 0)
                    {
                        process_tile(src_p, dst_p, word_idx - istart, rword_index - istart, length, 1);
                    }
                    else
                    {
                        process_tile(src_p, dst_p, word_idx - istart, rword_index - istart, length, -1);
                    }
                }
            }
        }

    private:

        static void read_tile(const SCA* __restrict data_ptr, SCA* __restrict tile_ptr, int stride)
        {
            for (size_t row=0; row < block_width; ++row)
            {
                int row_offset = row * stride;
                for (size_t col=0; col <block_width; ++col)
                {
                    *(tile_ptr++) = data_ptr[row_offset + col];
                }
            }
        }

        static void sign_tile(SCA tile[block_size], int sign) noexcept
        {
            for (size_t i = 0; i < block_size; ++i)
            {
                tile[i] = sign*tile[i]; // tile[i] = op(tile[i]);
            }
        }

        static void write_tile(SCA* __restrict tile_ptr, SCA* __restrict data_ptr, int stride)
        {
            for (size_t row=0; row<block_width; ++row)
            {
                int row_offset = row * stride;

                for (size_t col=0; col <block_width; ++col)
                {
                    data_ptr[row_offset + col] = *(tile_ptr++);
                }
            }
        }

        void process_tile(
                const SCA* input_data,
                SCA* output_data,
                size_t word_index,
                size_t rword_index,
                int degree,
                int sign
        ) const
        {
            SCA tile[block_size];
            auto stride = power(Width, degree+BlockLetters);

            read_tile(input_data + word_index*block_offset, tile, stride);

            reversing_permutation<Width, 2*BlockLetters> permutation;

            permutation(tile);

            sign_tile(tile, sign);

            write_tile(tile, output_data + rword_index*block_offset, stride);

        }
    };
} // dtl

template<typename Coeff>
class free_tensor_multiplication;

template<typename Coeff>
class shuffle_tensor_multiplication;

template<typename Coeff, DEG n_letters, DEG max_degree, typename...>
class shuffle_tensor;

template<typename Coeff>
class free_tensor_multiplication
{

    typedef typename Coeff::SCA scalar_t;

    template<typename Transform>
    class index_operator
    {
        Transform m_transform;

    public:
        index_operator(Transform t)
            : m_transform(t)
        {}

        void operator()(scalar_t* result_ptr, scalar_t const* lhs_ptr, scalar_t const* rhs_ptr, DIMN const lhs_target,
                        DIMN const rhs_target, bool assign = false)
        {
            scalar_t lhs;
            if (assign) {
                for (IDIMN i = 0; i < static_cast<IDIMN>(lhs_target); ++i) {
                    lhs = lhs_ptr[i];
                    for (IDIMN j = 0; j < static_cast<IDIMN>(rhs_target); ++j) {
                        *(result_ptr++) = m_transform(Coeff::mul(lhs, rhs_ptr[j]));
                    }
                }
            }
            else {
                for (IDIMN i = 0; i < static_cast<IDIMN>(lhs_target); ++i) {
                    lhs = lhs_ptr[i];
                    for (IDIMN j = 0; j < static_cast<IDIMN>(rhs_target); ++j) {
                        *(result_ptr++) += m_transform(Coeff::mul(lhs, rhs_ptr[j]));
                    }
                }
            }
        }
    };

    template<typename Transform>
    class key_operator
    {
        Transform m_transform;

    public:
        key_operator(Transform t)
            : m_transform(t)
        {}

        template<typename Vector>
        void
        operator()(Vector& result, typename Vector::KEY const& lhs_key, scalar_t const& lhs_val,
                   typename Vector::KEY const& rhs_key, scalar_t const& rhs_val)
        {
            result.add_scal_prod(lhs_key * rhs_key, m_transform(Coeff::mul(lhs_val, rhs_val)));
        }
    };

public:
    template<typename Algebra, typename Operator>
    Algebra& multiply_and_add(Algebra& result, Algebra const& lhs, Algebra const& rhs, Operator op) const
    {
        key_operator<Operator> kt(op);
        index_operator<Operator> it(op);
        lhs.buffered_apply_binary_transform(result, rhs, kt, it);
        return result;
    }

    template<typename Algebra, typename Operator>
    Algebra&
    multiply_and_add(Algebra& result, Algebra const& lhs, Algebra const& rhs, Operator op, DEG const max_depth) const
    {
        key_operator<Operator> kt(op);
        index_operator<Operator> it(op);
        lhs.buffered_apply_binary_transform(result, rhs, kt, it, max_depth);
        return result;
    }

    template<typename Algebra, typename Operator>
    Algebra multiply(Algebra const& lhs, Algebra const& rhs, Operator op) const
    {
        Algebra result;
        multiply_and_add(result, lhs, rhs, op);
        return result;
    }

    template<typename Algebra, typename Operator>
    Algebra multiply(Algebra const& lhs, Algebra const& rhs, Operator op, DEG const max_depth) const
    {
        Algebra result;
        multiply_and_add(result, lhs, rhs, op, max_depth);
        return result;
    }

    template<typename Algebra, typename Operator>
    Algebra& multiply_inplace(Algebra& lhs, Algebra const& rhs, Operator op) const
    {
        key_operator<Operator> kt(op);
        index_operator<Operator> it(op);
        lhs.unbuffered_apply_binary_transform(rhs, kt, it);
        return lhs;
    }

    template<typename Algebra, typename Operator>
    Algebra& multiply_inplace(Algebra& lhs, Algebra const& rhs, Operator op, DEG const max_depth) const
    {
        key_operator<Operator> kt(op);
        index_operator<Operator> it(op);
        lhs.unbuffered_apply_binary_transform(rhs, kt, it, max_depth);
        return lhs;
    }
};

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
                            free_tensor_basis<n_letters, max_degree>, Coeff, free_tensor_multiplication<Coeff>, VectorType, Args...>
{
    typedef free_tensor_multiplication<Coeff> multiplication_t;

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

    free_tensor& operator=(const free_tensor&) = default;
    //free_tensor& operator=(free_tensor&&) noexcept = default;

public:
    /// Ensures that the return type is a free_tensor.
    inline free_tensor operator*(const SCA& rhs) const
    {
        free_tensor result(*this);
        result *= rhs;
        return result;
    }

    /// Ensures that the return type is a free_tensor.
    inline free_tensor operator/(const RAT& rhs) const
    {
        free_tensor result(*this);
        result /= rhs;
        return result;
    }

    /// Ensures that the return type is a free_tensor.
    inline free_tensor operator*(const free_tensor& rhs) const
    {
        free_tensor result(*this);
        result *= rhs;
        return result;
    }

    /// Ensures that the return type is a free_tensor.
    inline free_tensor operator+(const free_tensor& rhs) const
    {
        free_tensor result(*this);
        result += rhs;
        return result;
    }

    /// Ensures that the return type is a free_tensor.
    inline free_tensor operator-(const free_tensor& rhs) const
    {
        free_tensor result(*this);
        result -= rhs;
        return result;
    }

    /// Ensures that the return type is a free_tensor.
    inline free_tensor operator-() const
    {
        return free_tensor(ALG::operator-());
    }

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

        for (auto cit = this->begin(); cit != this->end(); ++cit)
        {
            KEY temp_key = cit->key();
            KEY temp_key_reverse = temp_key.reverse();
            SCA temp_value = cit->value();

            int sign;

            if (temp_key.size() % 2 == 0)
            {
                sign = 1;
            }
            else
            {
                sign = -1;
            }

            result.add_scal_prod(temp_key_reverse, sign*temp_value);
        }

        return result;

    }

    // Implementation of the antipode for dense vector types.
    free_tensor antipode_impl(vectors::dtl::access_type_dense) const
    {
        free_tensor result;

#ifdef LIBALGEBRA_MAX_TILE_LETTERS
        constexpr DEG CalcLetters = integer_maths::logN(static_cast<unsigned>(LIBALGEBRA_L1_CACHE_SIZE),n_letters) / 2;
        constexpr DEG BlockLetters = (CalcLetters > LIBALGEBRA_MAX_TILE_LETTERS) ? LIBALGEBRA_MAX_TILE_LETTERS : CalcLetters;
#else
        constexpr DEG BlockLetters = integer_maths::logN(LIBALGEBRA_L1_CACHE_SIZE /sizeof(SCA), n_letters) / 2;
#endif

        const auto curr_degree = this->degree();

        vectors::dtl::vector_base_access::convert(result).resize_to_degree(curr_degree);

        // Get the pointers to the start of the data blob in memory.
        const SCA* src_ptr = vectors::dtl::data_access<VectorType<BASIS, Coeff>>::range_begin(vectors::dtl::vector_base_access::convert(*this));
        SCA* dst_ptr = vectors::dtl::data_access<VectorType<BASIS, Coeff>>::range_begin(vectors::dtl::vector_base_access::convert(result));

        dtl::tiled_inverse_operator<n_letters, max_degree, BlockLetters, SCA> t;

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

private:
    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive& ar, unsigned int const /* version */)
    {
        ar& boost::serialization::base_object<ALG>(*this);
    }
};

template<typename Coeff>
class shuffle_tensor_multiplication
{

    typedef typename Coeff::SCA scalar_t;

    /// Computes recursively the shuffle product of two keys
    template<typename Tensor>
    static Tensor _prod(typename Tensor::KEY const& k1, typename Tensor::KEY const& k2)
    {
        typedef typename Tensor::KEY key_type;

        typedef typename Tensor::BASIS basis_t;

        Tensor result;
        // unsigned i, j;
        const scalar_t one(+1);

        if ((basis_t::degree_tag::max_degree == 0) || (k1.size() + k2.size() <= basis_t::degree_tag::max_degree)) {
            if (k1.size() == 0) {
                result[k2] = one;
                return result;
            }
            if (k2.size() == 0) {
                result[k1] = one;
                return result;
            }
            // k1.size() >= 1 and k2.size() >= 1
            // let's just implement the multiplication
            const Tensor& first = prod<Tensor>(k1.rparent(), k2);
            const Tensor& second = prod<Tensor>(k1, k2.rparent());
            const key_type k1l{k1.lparent()}, k2l{k2.lparent()};

            typename Tensor::const_iterator cit;

            for (cit = first.begin(); cit != first.end(); ++cit) {
                result[k1l * cit->key()] += Coeff::one;
            }
            for (cit = second.begin(); cit != second.end(); ++cit) {
                result[k2l * cit->key()] += Coeff::one;
            }
        }

        return result;
    }

    /// The shuffle product of two basis elements.
    /**
    Returns the shuffle_tensor obtained by the concatenation product of two
    keys viewed as words of letters. The result is a unidimensional
    shuffle_tensor with a unique key (the concatenation of k1 and k2)
    associated to the +1 scalar. The already computed products are stored in
    a static mutiplication table to speed up further calculations.
    */
    template<typename Tensor>
    static const Tensor& prod(typename Tensor::KEY const& k1, typename Tensor::KEY const& k2)
    {
        typedef typename Tensor::KEY key_t;
        static boost::recursive_mutex table_access;
        // get exclusive recursive access for the thread
        boost::lock_guard<boost::recursive_mutex> lock(table_access);

        typedef std::map<std::pair<key_t, key_t>, Tensor> TABLE_T;
        static TABLE_T table;
        typename TABLE_T::iterator it;
        std::pair<key_t, key_t> p(std::min(k1, k2), std::max(k1, k2));
        it = table.find(p);
        if (it == table.end()) {
            return table[p] = _prod<Tensor>(k1, k2);
        }
        else {
            return it->second;
        }
    }

    template<typename Transform>
    class key_operator
    {
        Transform m_transform;

    public:
        explicit key_operator(Transform t)
            : m_transform(t)
        {}

        template<typename Vector>
        void
        operator()(Vector& result, typename Vector::KEY const& lhs_key, scalar_t const& lhs_val,
                   typename Vector::KEY const& rhs_key, scalar_t const& rhs_val)
        {
            result.add_scal_prod(prod<Vector>(lhs_key, rhs_key), m_transform(Coeff::mul(lhs_val, rhs_val)));
        }
    };

public:
    template<typename Algebra, typename Operator>
    Algebra& multiply_and_add(Algebra& result, Algebra const& lhs, Algebra const& rhs, Operator op) const
    {
        key_operator<Operator> kt(op);
        lhs.buffered_apply_binary_transform(result, rhs, kt);
        return result;
    }

    template<typename Algebra, typename Operator>
    Algebra&
    multiply_and_add(Algebra& result, Algebra const& lhs, Algebra const& rhs, Operator op, DEG const max_depth) const
    {
        key_operator<Operator> kt(op);
        lhs.buffered_apply_binary_transform(result, rhs, kt, max_depth);
        return result;
    }

    template<typename Algebra, typename Operator>
    Algebra multiply(Algebra const& lhs, Algebra const& rhs, Operator op) const
    {
        Algebra result;
        multiply_and_add(result, lhs, rhs, op);
        return result;
    }

    template<typename Algebra, typename Operator>
    Algebra multiply(Algebra const& lhs, Algebra const& rhs, Operator op, DEG const max_depth) const
    {
        Algebra result;
        multiply_and_add(result, lhs, rhs, op, max_depth);
        return result;
    }

    template<typename Algebra, typename Operator>
    Algebra& multiply_inplace(Algebra& lhs, Algebra const& rhs, Operator op) const
    {
        key_operator<Operator> kt(op);
        lhs.unbuffered_apply_binary_transform(rhs, kt);
        return lhs;
    }

    template<typename Algebra, typename Operator>
    Algebra& multiply_inplace(Algebra& lhs, Algebra const& rhs, Operator op, DEG const max_depth) const
    {
        key_operator<Operator> kt(op);
        lhs.unbuffered_apply_binary_transform(rhs, kt, max_depth);
        return lhs;
    }
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
                               shuffle_tensor_basis<n_letters, max_degree>, Coeff, shuffle_tensor_multiplication<Coeff>>
{
    typedef shuffle_tensor_multiplication<Coeff> multiplication_t;

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
    shuffle_tensor(LET
                           letter,
                   const SCA& s)
        : ALG(VECT::basis
                      .keyofletter(letter),
              s)
    {
    }

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

public:
    /// Ensures that the return type is a shuffle_tensor.
    inline shuffle_tensor operator*(const SCA& rhs) const
    {
        shuffle_tensor result(*this);
        result *= rhs;
        return result;
    }

    /// Ensures that the return type is a shuffle_tensor.
    inline shuffle_tensor operator/(const RAT& rhs) const
    {
        shuffle_tensor result(*this);
        result /= rhs;
        return result;
    }

    /// Ensures that the return type is a shuffle_tensor.
    inline shuffle_tensor operator*(const shuffle_tensor& rhs) const
    {
        shuffle_tensor result(*this);
        result *= rhs;
        return result;
    }

    /// Ensures that the return type is a shuffle_tensor.
    inline shuffle_tensor operator+(const shuffle_tensor& rhs) const
    {
        shuffle_tensor result(*this);
        result += rhs;
        return result;
    }

    /// Ensures that the return type is a shuffle_tensor.
    inline shuffle_tensor operator-(const shuffle_tensor& rhs) const
    {
        shuffle_tensor result(*this);
        result -= rhs;
        return result;
    }

    /// Ensures that the return type is a shuffle_tensor.
    inline shuffle_tensor operator-(void) const { return shuffle_tensor(ALG::operator-()); }

private:
    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive& ar, unsigned int const /* version */)
    {
        ar& boost::serialization::base_object<ALG>(*this);
    }
};

}// namespace alg
// Include once wrapper
#endif// DJC_COROPA_LIBALGEBRA_TENSORH_SEEN

// EOF.
