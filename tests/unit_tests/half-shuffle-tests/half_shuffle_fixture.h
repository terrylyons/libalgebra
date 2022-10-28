//
// Created by user on 08/06/22.
//

#ifndef LIBALGEBRA_TESTS_UNIT_TESTS_HALF_SHUFFLE_TESTS_HALF_SHUFFLE_FIXTURE_H_
#define LIBALGEBRA_TESTS_UNIT_TESTS_HALF_SHUFFLE_TESTS_HALF_SHUFFLE_FIXTURE_H_
#include <libalgebra/alg_types.h>
#include <libalgebra/alternative_multiplications.h>
#include <libalgebra/area_tensor_basis.h>
#include <libalgebra/area_tensor_multiplication.h>
#include <libalgebra/half_shuffle_tensor_basis.h>
#include <libalgebra/half_shuffle_tensor_multiplication.h>
#include <libalgebra/implementation_types.h>
#include <libalgebra/libalgebra.h>

#include <boost/range/algorithm.hpp>
#include <boost/range/iterator.hpp>
#include <boost/range/iterator_range.hpp>
#include <list>
#include <random>
#include <set>

template<size_t D, size_t W, coefficient_t F = Rational, vector_t VectorType = Hybrid>
struct HalfShuffleFixture : public alg_types<D, W, F, VectorType> {
    // the tensor and lie types
    typedef alg_types<D, W, Rational, VectorType> ALG_TYPES;

    using scalar_type = double;
    using integer_type = unsigned;

    // introduce a new class MDEGREE that can represent the mixed homogeneity type of a hall basis vector
    typedef typename alg::lie<typename ALG_TYPES::COEFF, ALG_TYPES::ALPHABET_SIZE, 1> MDEGREE;
    typedef typename ALG_TYPES::TENSOR TENSOR;
    typedef typename ALG_TYPES::SHUFFLE_TENSOR SHUFFLE_TENSOR;
    typedef typename alg::half_shuffle_tensor<typename ALG_TYPES::COEFF, ALG_TYPES::ALPHABET_SIZE, ALG_TYPES::DEPTH, alg::vectors::sparse_vector> HALF_SHUFFLE_TENSOR;
    typedef typename alg::area_tensor<typename ALG_TYPES::COEFF, ALG_TYPES::ALPHABET_SIZE, ALG_TYPES::DEPTH, alg::vectors::sparse_vector> AREA_TENSOR;
    typedef typename alg::free_tensor<typename ALG_TYPES::COEFF, ALG_TYPES::ALPHABET_SIZE, ALG_TYPES::DEPTH> FREE_TENSOR;

    //template <typename Tensor>
    //using free_multiply = alg::free_multiply<Tensor>;

    //template <typename Tensor>
    //using shuffle_multiply = alg::shuffle_multiply<Tensor>;

    //template <typename Tensor>
    //using half_shuffle_multiply =alg::half_shuffle_multiply<Tensor>;

    //template <typename Tensor>
    //using area_multiply = alg::area_multiply<Tensor>;

    typedef typename ALG_TYPES::MAPS MAPS;
    typedef typename ALG_TYPES::S S;
    typedef typename ALG_TYPES::DEG DEG;
    typedef typename ALG_TYPES::LIE LIE;
    typedef typename LIE::BASIS LBASIS;
    typedef typename LIE::KEY LKEY;
    typedef typename LBASIS::PARENT PARENTS;

    typedef typename TENSOR::BASIS BASIS;
    typedef typename TENSOR::KEY KEY;

    typedef typename SHUFFLE_TENSOR::BASIS SBASIS;
    typedef typename SHUFFLE_TENSOR::KEY SKEY;

    typedef typename HALF_SHUFFLE_TENSOR::BASIS HSBASIS;
    typedef typename HALF_SHUFFLE_TENSOR::KEY HSKEY;

    typedef typename AREA_TENSOR::BASIS ABASIS;
    typedef typename AREA_TENSOR::KEY AKEY;

    typedef typename alg::DIMN DIMN;
    //typedef typename alg::DEG DEG;

    // it is important that these are static
    // so they can be used in classes defined
    // within classes. However this syntax only works from 2017
    LBASIS& lbasis = LIE::basis;
    LBASIS& hall_set = lbasis;
    BASIS& basis = TENSOR::basis;
    SBASIS& sbasis = SHUFFLE_TENSOR::basis;
    HSBASIS& hsbasis = HALF_SHUFFLE_TENSOR::basis;

    //static bool letter(LKEY arg) noexcept
    //{
    //	return letter(arg);
    //}

    // the basic translation functions
    MAPS maps;
    // any function defined on letters {LKEY k | lbasis.letter(k) = true}
    // with values in a set with a binary operation extends to be defined
    // canonically on the full range of Hall set indexes 0 < k < hall_set.size().

    // the look up table for degree
    std::vector<S> degree;// the degree of each key (degree[k] should equal lbasis.degree(k))

    // the multi-degree of each key
    std::vector<MDEGREE> mdegree;

    // for each multi-degree the hall set indexes k that have that have that multi-degree
    std::multimap<MDEGREE, LKEY> typed_basis;

    // the tensor form of the lie basis element associated to each key
    std::vector<TENSOR> tensorlie;

    // the shuffle forms of the Hall areas associated to each key
    std::vector<SHUFFLE_TENSOR> hall_areas;
    ///////////////////////////
    ///// Constructs an instance from an algebra instance.
    //shuffle_tensor(const ALG& a) : ALG(a) {}

    ///// Constructs an instance from a sparse_vector instance.
    //shuffle_tensor(const VECT& v) : ALG(v) {}

    ///// Constructs a unidimensional instance from a letter and a scalar.
    //shuffle_tensor(LET
    //	letter,
    //	const SCA& s
    //)
    //	:
    //	ALG(VECT::basis
    //		.
    //		keyofletter(letter), s
    //	) {
    //}

    ///// Constructs a unidimensional instance from a key (basis element).
    //explicit shuffle_tensor(const KEY& k) : ALG(k) {}

    ///// Constructs a unidimensional instance from a scalar.
    //explicit shuffle_tensor(const SCA& s) : ALG(VECT::basis.empty_key, s) {}
    ///////////////////////////

    template<class OUT, class LEAF_INIT, class FUSE_PARENTS>
    void hall_fill(OUT out, LEAF_INIT leaf_init, FUSE_PARENTS fuse_parents) const
    {
        // degree should be empty
        out.clear();

        out.reserve(hall_set.size());

        // first entry in the Hall Set is not a hall element
        out.emplace_back();

        // populate the leaves directly and the remainder of the hall set indices via back reference.
        for (LKEY k1 = 1; k1 < hall_set.size(); ++k1)
            out.emplace_back(lbasis.letter(k1) ? leaf_init(k1) : fuse_parents(k1));
    }

    MDEGREE multi_degree(const LKEY k)
    {
        return mdegree[k];
    }

    static std::list<LKEY> UnPackWord(KEY word_key)
    {
        std::list<LKEY> word_lkey;
        KEY temp(word_key);
        while (temp.size()) {
            word_lkey.emplace_back(temp.FirstLetter());
            temp = temp.rparent();
        }
        return word_lkey;
    }

    std::list<LKEY> word_2_monotone_hall_seq(KEY word_key) const
    {
        // expand each word to a list of letters recorded as LKEY
        std::list<LKEY> word_lkey = UnPackWord(word_key);
        // rewrite the sequence until it is decreasing (in the hall sense - increasing as LKEYS!)

        if (!word_lkey.empty()) {
            typename std::list<LKEY>::iterator ptr = word_lkey.end(), ptr2;
            --ptr;
            while (ptr != word_lkey.begin()) {
                --(ptr2 = ptr);
                if (ptr != word_lkey.end() && *ptr2 > *ptr) {
                    // legal rise
                    // replace the two hall keys with one
                    PARENTS combined = std::make_pair(*ptr, *ptr2);
                    LKEY k;
                    try {
                        k = lbasis[combined];
                    }
                    catch (const std::exception& e) {
                        std::cout << " a standard exception '"
                                  << e.what() << "'\n";
                        std::cout << " the offending pair is "
                                  << combined.first << " " << combined.second << "\n";
                    }
                    auto iptr = word_lkey.erase(ptr2, ++ptr);
                    word_lkey.insert(iptr, k);
                }
                else
                    --ptr;
            }
        }
        return word_lkey;
    }

    class lk2f
    {
    public:
        lk2f(const LBASIS& lbasis)
            : liekey2foliage(lbasis, letter_to_foliage, liekey2foliage_binop)
        {
        }

    private:
        static KEY letter_to_foliage(const LKEY letter) noexcept
        {
            return KEY(letter);
        }
        /// THE HALL SET DEFINITION USED IN OUR CODE REVERSES THE TWO ENTRIES IN THE BRACKET
        /// So the foliage map needs to reverse them back.
        static KEY liekey2foliage_binop(const KEY a, const KEY b) noexcept
        {
            return b * a;
        }

        using liekey2foliage_type = typename LBASIS::template extended_function<
                decltype(&letter_to_foliage),
                decltype(&liekey2foliage_binop),
                //alg::lazy_cache_tag</*needs a predicATE TO INDICATE WHEN TO STOP*/>
                //alg::lookup_table_tag</*needs a predicATE TO INDICATE WHEN TO STOP*/>
                alg::no_caching_tag>;
        typedef typename liekey2foliage_type::output_type output_type;
        const liekey2foliage_type liekey2foliage;

    public:
        /// expand an lkey as the join of the letters in its parents
        typename liekey2foliage_type::output_type operator()(LKEY i) const
        {
            return liekey2foliage(i);
        }
    };

    class H2degree
    {
    public:
        H2degree(const LBASIS& lbasis)
            : hall2degree(lbasis, letter_to_degree, H2degree_binop)
        {
        }

    private:
        static DEG letter_to_degree(const LKEY letter) noexcept
        {
            return DEG(1);
        }

        static DEG H2degree_binop(const DEG a, const DEG b) noexcept
        {
            return a + b;
        }

        using H2degree_type = typename LBASIS::template extended_function<
                decltype(&letter_to_degree),
                decltype(&H2degree_binop),
                alg::lazy_cache_tag<void /*needs a predicATE TO INDICATE WHEN TO STOP*/>
                //alg::lookup_table_tag</*needs a predicATE TO INDICATE WHEN TO STOP*/>
                //alg::no_caching_tag
                >;
        typedef typename H2degree_type::output_type output_type;
        const H2degree_type hall2degree;

    public:
        /// expand an lkey as the join of the letters in its parents
        typename H2degree_type::output_type operator()(LKEY i) const
        {
            return hall2degree(i);
        }
    };

    class H2mdegree
    {
    public:
        H2mdegree(const LBASIS& lbasis)
            : hall2mdegree(lbasis, letter_to_mdegree, H2mdegree_binop)
        {
        }

    private:
        static MDEGREE letter_to_mdegree(const LKEY letter) noexcept
        {
            return MDEGREE(letter);
        }

        static MDEGREE H2mdegree_binop(const MDEGREE a, const MDEGREE b) noexcept
        {
            return a + b;
        }

        using H2mdegree_type = typename LBASIS::template extended_function<
                decltype(&letter_to_mdegree),
                decltype(&H2mdegree_binop),
                alg::lazy_cache_tag<void /*needs a predicATE TO INDICATE WHEN TO STOP*/>
                //alg::lookup_table_tag</*needs a predicATE TO INDICATE WHEN TO STOP*/>
                //alg::no_caching_tag
                >;
        typedef typename H2mdegree_type::output_type output_type;
        const H2mdegree_type hall2mdegree;

    public:
        /// expand an lkey as the join of the letters in its parents
        typename H2mdegree_type::output_type operator()(LKEY i) const
        {
            return hall2mdegree(i);
        }
    };

    class H2area
    {
    public:
        H2area(const LBASIS& hall_set)
            : hall2area(hall_set, letter_to_area, H2area_binop)
        {
        }

    private:
        static SHUFFLE_TENSOR letter_to_area(const SKEY sletter) noexcept
        {
            const S one(1);
            return SHUFFLE_TENSOR(sletter, one);
        }

        static SHUFFLE_TENSOR H2area_binop(const SHUFFLE_TENSOR& a, SHUFFLE_TENSOR& b) noexcept
        {
            return area_multiply(a, b);
        }

        using H2area_type = typename LBASIS::template extended_function<
                decltype(&letter_to_area),
                decltype(&H2area_binop),
                alg::lazy_cache_tag<void /*needs a predicATE TO INDICATE WHEN TO STOP*/>
                //alg::lookup_table_tag</*needs a predicATE TO INDICATE WHEN TO STOP*/ >
                //alg::no_caching_tag
                >;
        typedef typename H2area_type::output_type output_type;
        const H2area_type hall2area;

    public:
        /// expand an lkey or hall set element as the area of its parents
        typename H2area_type::output_type operator()(LKEY i) const
        {
            return hall2area(i);
        }
    };

    class H2half_shuffle
    {
    public:
        H2half_shuffle(const LBASIS& hall_set)
            : hall2half_shuffle(hall_set, letter_to_half_shuffle, H2half_shuffle_binop)
        {
        }

    private:
        static SHUFFLE_TENSOR letter_to_half_shuffle(const SKEY sletter) noexcept
        {
            return SHUFFLE_TENSOR(sletter, S(1));
        }

        static SHUFFLE_TENSOR H2half_shuffle_binop(const SHUFFLE_TENSOR& a, SHUFFLE_TENSOR& b) noexcept
        {
            return half_shuffle_multiply(a, b);
        }

        using H2half_shuffle_type = typename LBASIS::template extended_function<
                decltype(&letter_to_half_shuffle),
                decltype(&H2half_shuffle_binop),
                alg::lazy_cache_tag<void /*needs a predicATE TO INDICATE WHEN TO STOP*/>
                //alg::lookup_table_tag</*needs a predicATE TO INDICATE WHEN TO STOP*/ >
                //alg::no_caching_tag
                >;
        typedef typename H2half_shuffle_type::output_type output_type;
        const H2half_shuffle_type hall2half_shuffle;

    public:
        /// expand an lkey or hall set element as the half_shuffle of its parents
        typename H2half_shuffle_type::output_type operator()(LKEY i) const
        {
            return hall2half_shuffle(i);
        }
    };

    //class H2tensorlie
    //{
    //public:
    //	H2tensorlie(const LBASIS& lbasis, MAPS& maps_)
    //		:hall2tensorlie(lbasis, letter_to_tensorlie, H2tensorlie_binop) maps(maps_) {
    //	}
    //private:
    //	MAPS& maps;

    //	static TENSOR letter_to_tensorlie(const LKEY letter) noexcept {
    //		return maps.l2t(LIE(LKEY(letter)));
    //	}

    //	static TENSOR H2tensorlie_binop(const TENSOR& a, const TENSOR& b) noexcept {
    //		return a * b - b * a;
    //	}

    //	using H2tensorlie_type = typename LBASIS::template
    //		extended_function<
    //		decltype(&letter_to_tensorlie),
    //		decltype(&H2tensorlie_binop),
    //		//alg::lazy_cache_tag</*needs a predicATE TO INDICATE WHEN TO STOP*/>
    //		//alg::lookup_table_tag</*needs a predicATE TO INDICATE WHEN TO STOP*/>
    //		alg::no_caching_tag
    //		>;
    //	typedef typename H2tensorlie_type::output_type output_type;
    //	const H2tensorlie_type hall2tensorlie;

    //public:
    //	/// expand an lkey as the join of the letters in its parents
    //	typename H2tensorlie_type::output_type operator() (LKEY i) const {
    //		return hall2tensorlie(i);
    //	}
    //};

    HalfShuffleFixture()
        : liekey2foliage(lbasis)
    {

        auto degree_leaf_init = [](LKEY k) -> S { return S(1); };
        auto degree_fuse_parents = [&](LKEY k) -> S { return degree[hall_set[k].first] + degree[hall_set[k].second]; };
        hall_fill<typename std::remove_const<decltype(degree)>::type&, decltype(degree_leaf_init), decltype(degree_fuse_parents)>(degree, degree_leaf_init, degree_fuse_parents);

        auto mdegree_leaf_init = [](LKEY k) -> MDEGREE { return MDEGREE(k); };
        auto mdegree_fuse_parents = [&](LKEY k) -> MDEGREE { return mdegree[hall_set[k].first] + mdegree[hall_set[k].second]; };
        hall_fill<typename std::remove_const<decltype(mdegree)>::type&, decltype(mdegree_leaf_init), decltype(mdegree_fuse_parents)>(mdegree, mdegree_leaf_init, mdegree_fuse_parents);

        auto tensorlie_leaf_init = [&](LKEY k) -> TENSOR { return maps.l2t(LIE(LKEY(k))); };
        auto tensorlie_fuse_parents = [&](LKEY k) -> TENSOR { return tensorlie[hall_set[k].first] * tensorlie[hall_set[k].second]; };
        hall_fill<typename std::remove_const<decltype(tensorlie)>::type&, decltype(tensorlie_leaf_init), decltype(tensorlie_fuse_parents)>(tensorlie, tensorlie_leaf_init, tensorlie_fuse_parents);

        for (LKEY k = 1; k < hall_set.size(); ++k)
            typed_basis.emplace(mdegree[k], k);
    }

    // the signature kernel, coordinate wise
    S K(const TENSOR& lhs, const TENSOR& rhs) const
    {
        S ans = 0;
        for (auto it = begin(lhs); it != end(lhs); ++it)
            ans += it->value() * rhs[it->key()];
        return ans;
    }

    // decreasing hall index
    struct DHI : private std::vector<LKEY> {
        // index dependent data
        const MDEGREE _mdegree;
        const DEG _degree;

        // make the const data publicly available through the vector interface
        operator const std::vector<LKEY>&()
        {
            return *this;
        }

        // move constructor
        DHI(DHI&& other)
        noexcept
            : std::vector<LKEY>(std::move(static_cast<std::vector<LKEY>&>(other)))// is the move required?
              ,
              _mdegree(std::move(other._mdegree)), _degree(std::move(other._degree))
        // no guarantee that the degree and m degree are compatible with the sequence after the move
        {
        }

        // constructors from initialization lists of keys and vectors of keys
        template<class... T>
        DHI(T... args)
        noexcept
            : std::vector<LKEY>(args...)
        {
            std::sort(std::vector<LKEY>::begin(), std::vector<LKEY>::end());
            _mdegree = MDEGREE();
            for (auto k : static_cast<const std::vector<LKEY>&>(*this)) {
                _mdegree += mdegree[k];
            };
            _degree = _mdegree.NormL1();
        }

        DHI& operator=(const DHI& other) = delete;
    };

    void ShowHallSet() const
    {
        std::cout << "\nThe Hall set: " << std::endl;
        for (auto i = lbasis.begin(); i != lbasis.end(); i = lbasis.nextkey(i)) {
            std::cout << i << ": (" << lbasis[i].first << ", " << lbasis[i].second << ") ";
        }
        std::cout << std::endl;
    }

    std::list<LKEY>
    ToDecreasingIntegers(scalar_type test_word)
    {
        scalar_type param, fractpart, intpart;
        std::list<LKEY> expansion;
        param = test_word;
        while (true) {
            fractpart = modf(param, &intpart);
            if (fractpart <= 1. / (intpart + 1.))
                return expansion;
            else
                param = 1 / fractpart;
            expansion.push_back(LKEY(intpart));
        };
        return expansion;
    }

    template<typename Tensor>
    Tensor diff_shuffle_half_shuffle(const Tensor& lhs, const Tensor& rhs)
    {
        auto left = alg::shuffle_multiply(lhs, rhs);
        auto right = Tensor(lhs[KEY()]*rhs[KEY()]) + ((alg::half_shuffle_multiply(lhs, rhs))
                   + (alg::half_shuffle_multiply(rhs, lhs)));

        return left - right;
    }

    template<typename Tensor>
    Tensor diff_area_shuffle(const Tensor& lhs, const Tensor& rhs)
    {
        return ((alg::shuffle_multiply(lhs, rhs) + alg::area_multiply(lhs, rhs))
                - alg::half_shuffle_multiply(rhs, lhs) * S(2));
    }

private:
public:
    /// expand an lkey as the join of the letters in its parents
    const lk2f liekey2foliage;
};

#endif//LIBALGEBRA_TESTS_UNIT_TESTS_HALF_SHUFFLE_TESTS_HALF_SHUFFLE_FIXTURE_H_
