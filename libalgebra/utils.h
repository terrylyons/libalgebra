/* *************************************************************

Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurkó and Arend Janssen.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */

//  utils.h

// Include once wrapper
#ifndef DJC_COROPA_LIBALGEBRA_UTILSH_SEEN
#define DJC_COROPA_LIBALGEBRA_UTILSH_SEEN

#include <libalgebra/vector_bundle.h>

namespace alg {

/// Provides maps between lie and free_tensor instances.
template<typename Coeff, DEG n_letters, DEG max_degree, typename Tensor, typename Lie>
class maps
{
    typedef typename Coeff::S SCA;
    typedef typename Coeff::Q RAT;

    /// The Free Associative Algebra Basis type
    typedef free_tensor_basis<n_letters, max_degree> TBASIS;
    /// The Free Lie Associative Algebra Basis type
    typedef lie_basis<n_letters, max_degree> LBASIS;
    /// The Free Lie Associative Algebra Basis KEY type
    typedef typename LBASIS::KEY LKEY;
    /// The Free Associative Algebra Basis KEY type
    typedef typename TBASIS::KEY TKEY;
    /// The Free Lie Associative Algebra element type
    // typedef lie <SCA, RAT, n_letters, max_degree> LIE;
    /// The Free Associative Algebra element type
    // typedef free_tensor <SCA, RAT, n_letters, max_degree> TENSOR;
    typedef Lie LIE;
    typedef Tensor TENSOR;

    typedef vectors::sparse_vector<
            TBASIS, Coeff, std::unordered_map<TKEY, SCA, typename TKEY::hash>>
            sparse_tensor_vect;

    typedef algebra<
            TBASIS,
            Coeff,
            free_tensor_multiplication<n_letters, max_degree>,
            vectors::sparse_vector
    > sparse_tensor_t;

    using dense_lie1_t = algebra<
            lie_basis<n_letters, 1>,
            Coeff,
            lie_multiplication<n_letters, max_degree>,
            vectors::dense_vector
    >;

private:
    struct expand_letter {
        sparse_tensor_t operator()(const LET& l) const
        {
            return sparse_tensor_t(TKEY(l));
        }
    };

    struct commutator_type {
        sparse_tensor_t operator()(const sparse_tensor_t& a, const sparse_tensor_t& b) const
        {
            return commutator(a, b);
        }
    };

    using expand_function_t = typename LBASIS::template extended_function<expand_letter, commutator_type, lazy_cache_tag<void>>;

public:
    /**
     * @brief Returns the free_tensor corresponding to the Lie key k.
     * For performance reasons, the already computed expressions are stored in a
     * static table to speed up further calculus. The function returns a
     * constant reference to an element of this table.
     */
    expand_function_t expand;

    /// Default constructor.
    maps()
        : expand(LIE::basis.template extend_function<expand_function_t>())
    {}

public:
    /// computes the linear map
    class t2t
    {
        typedef alg::LET (*translator)(const LET);

        const translator h;

    public:
        explicit t2t(translator arg)
            : h(arg)
        {}

        template<typename Coeff2, DEG n_letters1, DEG max_degree1>
        TENSOR operator()(const alg::free_tensor<Coeff2, n_letters1, max_degree1>& in) const
        {
            typedef alg::free_tensor<Coeff2, n_letters1, max_degree1> TENSORIN;

            TENSOR out;
            for (typename TENSORIN::const_iterator it = in.begin(); it != in.end(); ++it) {
                typename TENSOR::KEY y(it->key(), h);
                if (SCA(0) == (out[y] += (it->value()))) {
                    out.erase(y);
                }
            }
            return out;
        }
    };

    /// Computes the free_tensor truncated exponential of a free lie element.
    inline TENSOR exp(LET l)
    {
        TKEY k;// empty word.
        TENSOR result(k);
        SCA coef(+1);
        DEG i;
        for (i = 1; i <= max_degree; ++i) {
            coef /= (RAT)i;
            k.push_back(l);
            result[k] = coef;
        }
        return result;
    }

    /// Returns the free_tensor corresponding to a free lie element.
    template<typename InputLie>
    inline Tensor l2t(const InputLie& arg)
    {
        Tensor result;
        typename InputLie::const_iterator i, iend(arg.end());
        for (i = arg.begin(); i != iend; ++i) {
            result.add_scal_prod(expand(i->key()), i->value());
        }
        return result;
    }

    template <typename InputLie>
    vector_bundle<Tensor> l2t(const vector_bundle<InputLie>& arg)
    {
        return {l2t(arg.base()), l2t(arg.fibre())};
    }

    /// Convert lie to tensor
    Tensor l2t(dense_lie1_t&& arg)
    {
        SCA* start = &arg.begin()->value();
        return Tensor(DIMN(1), start, start + n_letters);
    }

    /// Convert tensor to lie
    Tensor l2t(dense_lie1_t const& arg)
    {
        SCA const* start = &arg.begin()->value();
        return Tensor(DIMN(1), &*start, start + n_letters);
    }

    /**
     * @brief Returns the free lie element corresponding to a tensor_element.
    * This is the Dynkin map obtained by right bracketing. Of course, the
    * result makes sense only if the given free_tensor is the tensor expression
    * of some free lie element.
    */
    template<typename InputTensor>
    inline Lie t2l(const InputTensor& arg)
    {
        Lie result;
        typename InputTensor::const_iterator i;
        for (i = arg.begin(); i != arg.end(); ++i) {
            if (i->key().size() != 0 && i->value() != Tensor::zero) {
                result.add_scal_prod(rbraketing(i->key()), i->value());
            }
        }
        typename Lie::iterator j;
        for (j = result.begin(); j != result.end(); ++j) {
            j->value() /= (RAT)(LIE::basis.degree(j->key()));
        }
        return result;
    }

    template <typename InputTensor>
    vector_bundle<Lie> t2l(const vector_bundle<InputTensor>& arg)
    {
        return {t2l(arg.base()), t2l(arg.fibre())};
    }

    /**
     * @brief For a1,a2,...,an, return the expression [a1,[a2,[...,an]]].
     * For performance reasons, the already computed expressions are stored in a
     * static table to speed up further calculus. The function returns a
     * constant reference to an element of this table.
     */
    template<typename TensorKey>
    inline const LIE& rbraketing(const TensorKey& k)
    {
        //static boost::recursive_mutex table_access;
        // get exclusive recursive access for the thread
        //boost::lock_guard<boost::recursive_mutex> lock(table_access);

        //static boost::container::flat_map<TensorKey, LIE> lies;
        //typename boost::container::flat_map<TensorKey, LIE>::iterator it;
        static std::unordered_map<TensorKey, LIE, typename TensorKey::hash> lies;
        //static std::map<TensorKey, LIE> lies;
        //typename std::unordered_map<TensorKey, LIE, typename TensorKey::hash>::iterator it;
        LIE& value = lies[k];

        if (!value.empty()) {
            return value;
        }
        return value = _rbraketing(k);

        /*
        it = lies.find(k);
        if (it != lies.end()) {
            return it->second;
        } else {
            return lies[k] = _rbraketing(k);
        }
         */
    }

private:
    /// a1,a2,...,an is converted into [a1,[a2,[...,an]]] recursively.
    template<typename TensorKey>
    LIE _rbraketing(const TensorKey& k)
    {
        if (k.size() == 1) {
            return (LIE)LIE::basis.keyofletter(k.FirstLetter());
        }
        return rbraketing(k.lparent()) * rbraketing(k.rparent());
    }
};

/// Provides Campbell-Baker-Hausdorff formulas.
template<typename Coeff, DEG n_letters, DEG max_degree, typename Tensor, typename Lie>
class cbh
{
    typedef typename Coeff::S SCA;
    typedef typename Coeff::Q RAT;

    /// The Free Associative Algebra Basis type.
    typedef free_tensor_basis<n_letters, max_degree> TBASIS;
    /// The Free Lie Associative Algebra Basis type.
    typedef lie_basis<n_letters, max_degree> LBASIS;
    /// The Free Lie Associative Algebra Basis KEY type.
    typedef typename LBASIS::KEY LKEY;
    /// The Free Associative Algebra Basis KEY type.
    typedef typename TBASIS::KEY TKEY;
    /// The Free Lie Associative Algebra element type.
    // typedef lie <SCA, RAT, n_letters, max_degree> LIE;
    typedef Lie LIE;
    /// The Free Associative Algebra element type.
    // typedef free_tensor <SCA, RAT, n_letters, max_degree> TENSOR;
    typedef Tensor TENSOR;
    /// The MAPS type.
    typedef maps<Coeff, n_letters, max_degree, Tensor, Lie> MAPS;
    /// Maps between lie and free_tensor instances.
    mutable MAPS m_maps;// TJL added mutable
public:
    /// The empty free_tensor.
    TENSOR empty_tensor;
    /// The empty free lie element.
    LIE empty_lie;

public:
    /// Default constructor.
    cbh() = default;

public:
    /// Returns the CBH formula as a free lie element from a vector of letters.
    inline LIE basic(const std::vector<LET>& s) const
    {
        if (s.empty()) {
            return empty_lie;
        }
        TENSOR tmp(m_maps.exp(s[0]));
        typename std::string::size_type i;
        for (i = 1; i < s.size(); ++i) {
            tmp *= m_maps.exp(s[i]);
        }
        return m_maps.t2l(log(tmp));
    }

private:
    template <typename T>
    struct bundle_or_lie {
        using type = LIE;
    };

    template <typename T>
    struct bundle_or_lie<vector_bundle<T>> {
        using type = vector_bundle<LIE>;
    };

    template <typename T>
    using bundle_or_lie_t = typename bundle_or_lie<T>::type;

public:


//    /// Returns the CBH formula as a free lie element from an iterator to lie objects
//    template<typename InputIt>
////    std::enable_if_t<std::is_same<std::remove_cv_t<typename std::iterator_traits<InputIt>::value_type>, LIE>::value, LIE>
//    LIE
//    full(InputIt start, InputIt finish)
//    {
//        if (start == finish) {
//            return empty_lie;
//        }
//
//        InputIt it(start);
//        TENSOR result(exp(m_maps.l2t(*(it++))));
//
//        for (; it != finish; ++it) {
//            result.fmexp_inplace(m_maps.l2t(*it));
//        }
//
//        return m_maps.t2l(log(result));
//    }

    /// Returns the CBH formula as a free lie element from an iterator to lie objects
    template<typename InputIt>
//    std::enable_if_t<std::is_same<std::remove_cv_t<typename std::iterator_traits<InputIt>::value_type>, vector_bundle<LIE>>::value, vector_bundle<LIE>>
    bundle_or_lie_t<typename std::iterator_traits<InputIt>::value_type>
    full(InputIt start, InputIt finish)
    {
        if (start == finish) {
            return {};
        }

        InputIt it(start);
        auto result = exp(m_maps.l2t(*(it++)));

        for (; it != finish; ++it) {
            result.fmexp_inplace(m_maps.l2t(*it));
        }

        return m_maps.t2l(log(result));
    }




    /// Returns the CBH formula as a free lie element from a vector of lie.
    inline LIE full(const std::vector<const LIE*>& lies) const
    {
        if (lies.empty()) {
            return empty_lie;
        }
        typename std::vector<const LIE*>::size_type i;
        TENSOR tmp(exp(m_maps.l2t(*lies[0])));
        for (i = 1; i < lies.size(); ++i) {
            tmp.fmexp_inplace(m_maps.l2t(*lies[i]));
        }
        return m_maps.t2l(log(tmp));
    }


    vector_bundle<Lie> full(const std::vector<const vector_bundle<Lie>*> lies) const
    {
        if (lies.empty()) {
            return vector_bundle<Lie>();
        }

        auto acc = exp(m_maps.l2t(*lies[0]));
        for (DIMN i=1; i<lies.size(); ++i) {
            acc.fmexp_inplace(m_maps.l2t(*lies[i]));
        }

        return m_maps.t2l(log(acc));
    }

};

}// namespace alg
// Include once wrapper
#endif// DJC_COROPA_LIBALGEBRA_UTILSH_SEEN

// EOF.
