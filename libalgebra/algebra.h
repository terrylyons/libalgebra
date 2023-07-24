/* *************************************************************

Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurkó and Arend Janssen.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */

//  algebra.h

// Include once wrapper
#ifndef DJC_COROPA_LIBALGEBRA_ALGEBRAH_SEEN
#define DJC_COROPA_LIBALGEBRA_ALGEBRAH_SEEN

#include "implementation_types.h"

#include <algorithm>
#include <functional>
#include <map>
#include <mutex>
#include <type_traits>
#include <unordered_map>
#include <utility>

#include <boost/call_traits.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/functional/hash.hpp>
#include <boost/type_traits/type_identity.hpp>


#include "basis.h"
#include "detail/meta.h"
#include "hybrid_vector.h"
#include "multiplication_helpers.h"
#include "tags.h"
#include "vectors.h"


namespace alg {

namespace dtl {

template<typename Tag>
struct has_degree_impl : public std::false_type {
};

template<DEG D>
struct has_degree_impl<basis::with_degree<D>>
    : public std::enable_if_t<(D > 0), std::true_type> {
};

template<typename Vector>
using vector_deg_tag = typename basis::basis_traits<
        typename Vector::BASIS::degree_tag>::degree_tag;

template<typename Vector>
using has_degree_t = has_degree_impl<vector_deg_tag<Vector>>;

template <typename LhsType, typename RhsType, typename OutType=LhsType>
struct original_algebras
{
    using lhs_type = LhsType;
    using rhs_type = RhsType;
    using out_type = OutType;

    const LhsType& lhs;
    const RhsType& rhs;
    OutType& out;
};


template<typename Multiplication>
class multiplication_traits
{
    using mul_reference = const Multiplication&;

    template<typename Result, typename LeftVector, typename RightVector, typename Fn>
    static void fma_impl(mul_reference mul, Result& out, const LeftVector& lhs, const RightVector& rhs, Fn op)
    {
        original_algebras<LeftVector, RightVector, Result> orig {lhs, rhs, out};
        mul.fma(out.base_vector(), lhs.base_vector(), rhs.base_vector(), op, orig);
    }

    template<typename Result, typename LeftVector, typename RightVector, typename Fn>
    static void fma_impl(mul_reference mul, Result& out, const LeftVector& lhs, const RightVector& rhs, Fn op, DEG max_degree)
    {
        original_algebras<LeftVector, RightVector, Result> orig {lhs, rhs, out};
        mul.fma(out.base_vector(), lhs.base_vector(), rhs.base_vector(), op, max_degree, orig);
    }

    template<typename LeftVector, typename RightVector, typename Fn>
    static void multiply_inplace_impl(mul_reference mul, LeftVector& lhs, const RightVector& rhs, Fn op)
    {
        original_algebras <LeftVector, RightVector> orig {lhs, rhs, lhs};
        mul.multiply_inplace(lhs.base_vector(), rhs.base_vector(), op, orig);
    }

    template<typename LeftVector, typename RightVector, typename Fn>
    static void multiply_inplace_impl(mul_reference mul, LeftVector& lhs, const RightVector& rhs, Fn op, DEG max_degree)
    {
        original_algebras<LeftVector, RightVector> orig {lhs, rhs, lhs};
        mul.multiply_inplace(lhs.base_vector(), rhs.base_vector(), op, max_degree, orig);
    }

    template<typename Result, typename LeftVector, typename RightVector, typename Fn>
    static void dispatch(mul_reference mul, Result& out, const LeftVector& lhs, const RightVector& rhs,
                         Fn op, basis::without_degree)
    {
        fma_impl(mul, out, lhs, rhs, op);
    }

    template<typename Result, typename LeftVector, typename RightVector, typename Fn, DEG D>
    static void dispatch(mul_reference mul, Result& out, const LeftVector& lhs, const RightVector& rhs,
                         Fn op, basis::with_degree<D>)
    {
        fma_impl(mul, out, lhs, rhs, op, D);
    }

    template<typename LeftVector, typename RightVector, typename Fn>
    static void dispatch(mul_reference mul, LeftVector& lhs, const RightVector& rhs,
                         Fn op, basis::without_degree)
    {
        multiply_inplace_impl(mul, lhs, rhs, op);
    }

    template<typename LeftVector, typename RightVector, typename Fn, DEG D>
    static void dispatch(mul_reference mul, LeftVector& lhs, const RightVector& rhs,
                         Fn op, basis::with_degree<D>)
    {
        multiply_inplace_impl(mul, lhs, rhs, op, D);
    }

    template<typename Result, typename LeftVector, typename RightVector, typename Fn>
    static void dispatch(mul_reference mul, Result& out, const LeftVector& lhs,
                         const RightVector& rhs, Fn op)
    {
        dispatch(mul, out, lhs, rhs, op, Result::degree_tag);
    }

    template<typename LeftVector, typename RightVector, typename Fn>
    static void dispatch(mul_reference mul, LeftVector& lhs, const RightVector& rhs,
                         Fn op)
    {
        dispatch(mul, lhs, rhs, op, LeftVector::degree_tag);
    }

public:
    template<typename Result,
             typename LeftVector,
             typename RightVector,
             typename Fn>
    static void multiply_and_add(mul_reference mul,
                                 Result& result,
                                 const LeftVector& lhs,
                                 const RightVector& rhs,
                                 Fn op)
    {
        dispatch(mul, result, lhs, rhs, op);
    }

    template<typename Result,
             typename LeftVector,
             typename RightVector>
    static void multiply_and_add(mul_reference mul,
                                 Result& result,
                                 const LeftVector& lhs,
                                 const RightVector& rhs)
    {
        dispatch(mul, result, lhs, rhs, mult::scalar_passthrough());
    }

    template<typename LeftVector, typename RightVector, typename Fn>
    static void multiply_inplace(mul_reference mul, LeftVector& lhs, const RightVector& rhs, Fn op)
    {
        dispatch(mul, lhs, rhs, op);
    }

    template<typename LeftVector, typename RightVector>
    static void multiply_inplace(mul_reference mul, LeftVector& lhs, const RightVector& rhs)
    {
        dispatch(mul, lhs, rhs, mult::scalar_passthrough());
    }

    template<typename Result,
             typename LeftVector,
             typename RightVector,
             typename Fn>
    static void multiply_and_add(mul_reference mul,
                                 Result& result,
                                 const LeftVector& lhs,
                                 const RightVector& rhs,
                                 Fn op,
                                 DEG max_degree)
    {
        fma_impl(mul, result, lhs, rhs, op, max_degree);
    }

    template<typename LeftVector, typename RightVector, typename Fn>
    static void multiply_inplace(mul_reference mul, LeftVector& lhs,
                                 const RightVector& rhs, Fn op,
                                 DEG max_degree)
    {
        multiply_inplace_impl(mul, lhs, rhs, op, max_degree);
    }
};

template<typename Vector>
class multiplication_helper
{
public:
    using basis_type = typename Vector::BASIS;
    using key_type = typename basis_type::KEY;
    using coeff_ring = typename Vector::coefficient_ring;
    using scalar_type = typename coeff_ring::S;
    using key_value = std::pair<key_type, scalar_type>;

protected:
    std::vector<key_value> buffer;

public:
    using const_iterator = typename std::vector<key_value>::const_iterator;

    explicit multiplication_helper(const Vector& rhs)
    {
        buffer.reserve(rhs.size());
        for (auto item : rhs) {
            buffer.push_back({item.key(), item.value()});
        }
    }

    const_iterator begin() const noexcept { return buffer.begin(); }
    const_iterator end() const noexcept { return buffer.end(); }
};

template<typename It>
class degree_range_iterator
{
    It m_begin, m_end;

public:
    using iterator = It;

    degree_range_iterator(It begin, It end)
        : m_begin(begin), m_end(end)
    {}

    iterator begin() noexcept { return m_begin; }
    iterator end() noexcept { return m_end; }
};

template<typename Vector>
class graded_multiplication_helper : public multiplication_helper<Vector>
{
    using Self = graded_multiplication_helper;
    using base_iter = typename std::vector<typename Self::key_value>::const_iterator;
    using base = multiplication_helper<Vector>;

protected:
    std::vector<base_iter> degree_ranges;
    DEG m_depth;

public:
    using range_iterator = degree_range_iterator<base_iter>;

    graded_multiplication_helper(const Vector& rhs, DEG max_degree)
        : multiplication_helper<Vector>(rhs)
    {
        using traits = basis::basis_traits<typename Self::basis_type>;
        typename traits::ordering_tag::pair_order ordering;

        std::sort(Self::buffer.begin(), Self::buffer.end(), ordering);
        const auto& basis = rhs.basis;

        m_depth = std::min(max_degree, traits::degree_tag::max_degree);
        auto it = base::begin();
        auto end = base::end();

        degree_ranges.resize(m_depth + 1, base::end());
        DEG degree = 0;
        for (; it != end; ++it) {
            auto current = basis.degree(it->first);
            while (degree < current) {
                degree_ranges[degree++] = it;
            }
        }

//        assert(degree_ranges.back() == base::buffer.end());
    }

    DEG depth() const noexcept { return m_depth; }

    range_iterator degree_range(DEG degree) const noexcept
    {
        if (degree > m_depth) {
            return {base::end(), base::end()};
        }
        return {base::begin(), degree_ranges[degree]};
    }
};

}// namespace dtl

/**
 * @brief Base class for multilpiers
 * @tparam Multiplier
 * @tparam Basis
 * @tparam Scalar
 * @tparam SSO
 */
template<typename Multiplier, typename Basis, DEG SSO = 1, typename Scalar = int>
class multiplier_base
{
protected:
    using basis_type = Basis;
    using key_type = typename Basis::KEY;
    using pair_type = std::pair<key_type, Scalar>;
    using inner_result_type = boost::container::small_vector<pair_type, SSO>;
    using result_type = const boost::container::small_vector_base<pair_type>&;
    using argument_type = typename  boost::call_traits<key_type>::param_type;

    struct cache_entry {
        inner_result_type data;
        bool set = false;
    };
    using cache_type = std::unordered_map<std::pair<key_type, key_type>,
                                          cache_entry,
                                          boost::hash<std::pair<key_type, key_type>>>;

    static inner_result_type fill(const std::map<key_type, Scalar>& arg)
    {
        inner_result_type result;
        result.reserve(arg.size());
        for (const auto& item : arg) {
            result.emplace_back(item.first, item.second);
        }
        return result;
    }

    static inner_result_type uminus(result_type arg)
    {
        inner_result_type result;
        result.reserve(arg.size());
        for (const auto& item : arg) {
            result.emplace_back(item.first, -item.second);
        }
        return result;
    }
    static inner_result_type add(result_type lhs, result_type rhs)
    {
        std::map<key_type, Scalar> tmp(lhs.begin(), lhs.end());

        for (const auto& item : rhs) {
            tmp[item.first] += item.second;
        }
        return fill(tmp);
        //        return {tmp.begin(), tmp.end()};
    }
    static inner_result_type sub(result_type lhs, result_type rhs)
    {
        std::map<key_type, Scalar> tmp(lhs.begin(), lhs.end());
        for (const auto& item : rhs) {
            tmp[item.first] -= item.second;
        }
        return fill(tmp);
        //        return {tmp.begin(), tmp.end()};
    }
    inner_result_type mul(key_type lhs, result_type rhs) const
    {
        std::map<key_type, Scalar> tmp;

        const auto& mul = static_cast<const Multiplier&>(*this);

        for (const auto& outer : rhs) {
            for (const auto& inner : mul(lhs, outer.first)) {
                tmp[inner.first] += inner.second * outer.second;
            }
        }
        return fill(tmp);
        //        return {tmp.begin(), tmp.end()};
    }
    inner_result_type mul(result_type lhs, key_type rhs) const
    {
        std::map<key_type, Scalar> tmp;
        const auto& mul = static_cast<const Multiplier&>(*this);
        for (const auto& outer : lhs) {
            for (const auto& inner : mul(outer.first, rhs)) {
                tmp[inner.first] += outer.second * inner.second;
            }
        }
        return fill(tmp);
    }

    result_type cached_compute(argument_type lhs, argument_type rhs) const
    {
        static std::recursive_mutex lock;
        static cache_type cache;

        std::lock_guard<std::recursive_mutex> access(lock);

        auto& result = cache[std::make_pair(lhs, rhs)];
        if (result.set) {
            return result.data;
        }
        result.set = true;
        return result.data = static_cast<const Multiplier&>(*this).prod_impl(lhs, rhs);
    }
};

namespace dtl {

template<typename Derived>
struct derived_or_this {
    template<typename Base>
    static const Derived& cast(const Base& arg) noexcept
    {
        return static_cast<const Derived&>(arg);
    }
};

template<>
struct derived_or_this<void> {
    template<typename Base>
    static const Base& cast(const Base& arg) noexcept
    {
        return arg;
    }
};

}// namespace dtl

/*
 * The new multiplication structure is as follows: Every multiplication is given
 * a unique class that provides methods for both inplace and external fused multiply
 * add operations. These operations should be templated over all vector arguments,
 * and specialised where optimisations exist (e.g. free tensor multiplication).
 *
 * The general case is handled by a multiplier, which is a function object
 * (with operator()) taking two key_type arguments and returning something akin to
 * a vector of key-value pairs that can be incorporated into the result vector.
 *
 *
 * The multiplier object is then used to declare a multiplication object described
 * above, with the basic implementations coming from base_multiplication below.
 */
template<typename Multiplier, typename Derived = void>
class base_multiplication
{

    using caster = dtl::derived_or_this<Derived>;

protected:
    Multiplier m_multiplier;

public:
    using multiplier_type = Multiplier;
    using basis_type = typename Multiplier::basis_type;

private:
    template<typename T>
    using is_compatible = typename vectors::vector_traits<T>::template is_compatible<basis_type>;

    template<typename V1, typename V2, typename V3 = void>
    using checked_out_t = std::enable_if_t<
            is_compatible<V1>::value && is_compatible<V2>::value && utils::void_or<V3, is_compatible, std::true_type>::type::value>;

public:
    template<typename OutVector, typename Add, typename Scalar>
    void asp_helper(OutVector& out, Add&& add, Scalar scal) const
    {
        using scalar_type = typename OutVector::SCALAR;
        for (const auto& item : add) {
            //            std::cout << std::make_pair(&out.basis, item.first) << ' ' << item.second << '\n';
            out.add_scal_prod(item.first, scalar_type(item.second) * scalar_type(scal));
        }
    }

    typename Multiplier::result_type
    eval(typename Multiplier::argument_type lhs, typename Multiplier::argument_type rhs) const
    {
        return m_multiplier(lhs, rhs);
    }

    template<typename OutVector, typename LeftVector, typename RightVector, typename Fn, typename O, typename=std::enable_if_t<!std::is_integral<O>::value>>
    checked_out_t<OutVector, LeftVector, RightVector>
    fma(OutVector& out, const LeftVector& left, const RightVector& right, Fn op, O) const
    {
        dtl::multiplication_helper<RightVector> helper(right);

        for (auto litem : left) {
            for (const auto& ritem : helper) {
                asp_helper(out, m_multiplier(litem.key(), ritem.first),
                           op(litem.value() * ritem.second));
            }
        }
    }

    template<typename OutVector, typename LeftVector, typename RightVector, typename Fn, typename O>
    checked_out_t<OutVector, LeftVector, RightVector>
    fma(OutVector& out, const LeftVector& left, const RightVector& right, Fn op, DEG max_degree, O) const
    {
        using out_traits = basis::basis_traits<typename OutVector::BASIS>;
        const auto max_d = std::min(max_degree, out_traits::degree_tag::max_degree);

        dtl::graded_multiplication_helper<RightVector> helper(right, max_d);

        const auto out_deg = std::min(max_d, left.degree() + right.degree());
        const auto& basis = out.basis;

        for (auto litem : left) {
            auto lkey = litem.key();
            auto lhs_deg = basis.degree(lkey);
            auto rhs_deg = out_deg - lhs_deg;
            for (const auto& ritem : helper.degree_range(rhs_deg)) {
                asp_helper(out, m_multiplier(lkey, ritem.first),
                           op(litem.value() * ritem.second));
            }
        }
    }

    template<typename LeftVector, typename RightVector, typename Fn, typename O, typename = std::enable_if_t<!std::is_integral<O>::value>>
    checked_out_t<LeftVector, RightVector> multiply_inplace(LeftVector& left, const RightVector& right, Fn op, O orig) const
    {
        if (!right.empty() && !left.empty()) {
            LeftVector tmp;
            caster::cast(*this).fma(tmp, left, right, op, orig);
            left.swap(tmp);
        }
        else {
            left.clear();
        }
    }

    template<typename LeftVector, typename RightVector, typename Fn, typename O>
    checked_out_t<LeftVector, RightVector> multiply_inplace(LeftVector& left, const RightVector& right, Fn op, DEG max_degree, O orig) const
    {
        if (!right.empty() && !left.empty()) {
            LeftVector tmp;
            caster::cast(*this).fma(tmp, left, right, op, max_degree, orig);
            left.swap(tmp);
        }
        else {
            left.clear();
        }
    }
};

namespace dtl {

/**
 * @brief Mixin to provide optimisation for the hybrid vector (for now)
 * @tparam Multiplication
 */
template<typename Multiplication>
class hybrid_vector_mixin
{

    template<typename Basis, typename Coeffs, typename Fn>
    void fma_mixed(
            vectors::hybrid_vector<Basis, Coeffs>& out,
            const vectors::hybrid_vector<Basis, Coeffs>& lhs,
            const vectors::hybrid_vector<Basis, Coeffs>& rhs,
            Fn op) const
    {
        using sparse_vec = vectors::sparse_vector<Basis, Coeffs>;

        if (lhs.sparse_empty() && rhs.sparse_empty()) {
            return;
        }

        const auto& self = static_cast<const Multiplication&>(*this);

        multiplication_helper<sparse_vec> helper(rhs.sparse_part());

        // dense * sparse
        if (lhs.dense_dimension() > 0) {
            auto lbegin = lhs.dense_part().begin();
            auto lend = lhs.dense_part().end();

            for (auto lit = lbegin; lit != lend; ++lit) {
                for (const auto& ritem : helper) {
                    self.asp_helper(out, self.eval(lit->key(), ritem.first),
                                    op(lit->value() * ritem.second));
                }
            }
        }

        // sparse * dense and sparse * sparse together
        auto lit = lhs.sparse_begin();
        auto lend = lhs.sparse_end();

        auto rbegin = rhs.dense_part().begin();
        auto rend = rhs.dense_part().end();
        auto rhs_has_dense = rhs.dense_dimension() > 0;

        for (; lit != lend; ++lit) {
            if (rhs_has_dense) {
                for (auto rit = rbegin; rit != rend; ++rit) {
                    self.asp_helper(out, self.eval(lit->key(), rit->key()),
                                    op(lit->value() * rit->value()));
                }
            }

            for (const auto& item : helper) {
                self.asp_helper(out, self.eval(lit->key(), item.first),
                                op(lit->value() * item.second));
            }
        }
        out.maybe_resize();
    }

    template<typename Basis, typename Coeffs, typename Fn>
    void fma_mixed(
            vectors::hybrid_vector<Basis, Coeffs>& out,
            const vectors::hybrid_vector<Basis, Coeffs>& lhs,
            const vectors::hybrid_vector<Basis, Coeffs>& rhs,
            Fn op,
            DEG max_degree) const
    {
        using sparse_vec = vectors::sparse_vector<Basis, Coeffs>;
        const auto& basis = out.basis;
        const auto& self = static_cast<const Multiplication&>(*this);

        if (lhs.sparse_empty() && rhs.sparse_empty()) {
            return;
        }

        graded_multiplication_helper<sparse_vec> helper(rhs, max_degree);
        auto out_deg = std::min(helper.depth(), lhs.degree() + rhs.degree());

        // dense*sparse
        if (lhs.dense_dimension() > 0) {
            auto ld_degree = lhs.dense_degree();
            auto rbegin = helper.begin();

            for (DEG ldeg = 0; ldeg <= ld_degree; ++ldeg) {
                auto rend = helper.degree_range(out_deg - ldeg).end();
                for (auto idx = basis.start_of_degree(ldeg); idx < basis.start_of_degree(ldeg + 1); ++idx) {
                    for (auto rit = rbegin; rit != rend; ++rit) {
                        self.asp_helper(out, self.eval(basis.index_to_key(idx), rit->first),
                                        op(lhs.dense_value(idx) * rit->second));
                    }
                }
            }
        }

        // sparse*dense and sparse*sparse
        auto lbegin = lhs.sparse_part().begin();
        auto lend = lhs.sparse_part().end();
        auto rddeg = rhs.dense_degree();
        auto rhs_has_dense = rhs.dense_dimension() > 0;

        for (auto lit = lbegin; lit != lend; ++lit) {
            auto lkey = lit->key();
            auto lhs_deg = basis.degree(lkey);
            assert(out_deg >= lhs_deg);
            auto rhs_deg = out_deg - lhs_deg;
            auto rhs_ddeg = std::min(rhs_deg, rddeg);

            if (rhs_has_dense) {
                for (DIMN idx = 0; idx < basis.start_of_degree(rhs_ddeg + 1); ++idx) {
                    self.asp_helper(out, self.eval(lkey, basis.index_to_key(idx)), op(lit->value() * rhs.dense_value(idx)));
                }
            }

            for (const auto& ritem : helper.degree_range(out_deg - lhs_deg)) {
                self.asp_helper(out, self.eval(lkey, ritem.first), op(lit->value() * ritem.second));
            }
        }

        out.maybe_resize();
    }

public:
    template<typename Basis, typename Coeffs, typename Fn>
    void fma(vectors::hybrid_vector<Basis, Coeffs>& out,
             const vectors::hybrid_vector<Basis, Coeffs>& lhs,
             const vectors::hybrid_vector<Basis, Coeffs>& rhs,
             Fn op) const
    {
        const auto& self = static_cast<const Multiplication&>(*this);

        self.fma(out.dense_part(), lhs.dense_part(), rhs.dense_part(), op);
        fma_mixed(out, lhs, rhs, op);
    }

    template<typename Basis, typename Coeffs, typename Fn>
    void fma(vectors::hybrid_vector<Basis, Coeffs>& out,
             const vectors::hybrid_vector<Basis, Coeffs>& lhs,
             const vectors::hybrid_vector<Basis, Coeffs>& rhs,
             Fn op,
             DEG max_degree) const
    {
        const auto& self = static_cast<const Multiplication&>(*this);
        self.fma(out.dense_part(), lhs.dense_part(), rhs.dense_part(), op, max_degree);
        //        std::cout << "BEFORE " << out << '\n';
        fma_mixed(out, lhs, rhs, op, max_degree);
        //        std::cout << "AFTER " << out << '\n';
    }

    template<typename Basis, typename Coeffs, typename Fn>
    void multiply_inplace(vectors::hybrid_vector<Basis, Coeffs>& lhs, const vectors::hybrid_vector<Basis, Coeffs>& rhs, Fn op) const
    {
        vectors::hybrid_vector<Basis, Coeffs> tmp;
        fma(tmp, lhs, rhs, op);
        lhs.swap(tmp);
    }

    template<typename Basis, typename Coeffs, typename Fn>
    void multiply_inplace(vectors::hybrid_vector<Basis, Coeffs>& lhs, const vectors::hybrid_vector<Basis, Coeffs>& rhs, Fn op, DEG max_degree) const
    {
        vectors::hybrid_vector<Basis, Coeffs> tmp;
        fma(tmp, lhs, rhs, op, max_degree);
        lhs.swap(tmp);
    }
};




}// namespace dtl






/**
 * @brief Class to store and manipulate associative algebra elements
 *
 * An algebra is a vector space that is given a compatible multiplication operation.
 * An algebra class is constructed in the same way. We provide the basis and coefficients
 * required to instantiate a vector class and provide an additional multiplication class
 * that implements the multiplication via key operators (acting on sparse data) and index
 * operators (operating on dense data).
 *
 * @tparam Basis Basis of underlying vector space
 * @tparam Coeff Coefficient field of underlying vector space
 * @tparam Multiplication Multiplication operation
 * @tparam VectorType Underlying vector data type to use; e.g. dense_vector or sparse_vector
 */
template<typename Basis, typename Coeff, typename Multiplication,
         template<typename, typename, typename...> class VectorType = alg::vectors::template_vector_type_selector<Basis, Coeff>::template type,
         typename Derived = void,
         typename... Args>
class algebra : public vectors::vector<Basis, Coeff, VectorType, Args...>
{

    typedef mult::scalar_passthrough scalar_passthrough;
    typedef mult::scalar_minus<Coeff> scalar_minus;
    typedef mult::scalar_post_mult<Coeff> scalar_post_mult;
    typedef mult::rational_post_div<Coeff> rational_post_div;

//    using derived_type = typename utils::void_or<Derived, boost::type_identity, algebra>::type;

    using derived_type = std::conditional_t<std::is_void<Derived>::value, algebra, Derived>;

public:
    typedef Basis BASIS;
    /// The inherited sparse vector type.
    typedef vectors::vector<Basis, Coeff, VectorType, Args...> VECT;
    /// Import of the iterator type from sparse_vector.
    typedef typename VECT::iterator iterator;
    /// Import of the constant iterator type from sparse_vector.
    typedef typename VECT::const_iterator const_iterator;
    /// Import of the KEY type from sparse_vector.
    typedef typename VECT::KEY KEY;
    /// Import of the SCALAR type from sparse_vector.
    typedef typename VECT::SCALAR SCALAR;
    /// Import of the RATIONAL type from sparse_vector.
    typedef typename VECT::RATIONAL RATIONAL;
    using typename VECT::coefficient_field;

    typedef Multiplication multiplication_t;

private:
    static const multiplication_t s_multiplication;
    using mtraits = dtl::multiplication_traits<Multiplication>;

public:
    /// Default constructor.
    /**
    Constructs an empty algebra element.
    */
    algebra()
        : VECT()
    {}

    /// Copy constructor.
    algebra(const algebra& a)
        : VECT(a)
    {}

    /// Move constructor
    algebra(algebra&& a) noexcept
        : VECT(std::move(a))
    {}

    /// Constructs an algebra instance from a sparse_vector.
    explicit algebra(const VECT& v)
        : VECT(v)
    {}

    /// Unidimensional constructor.
    explicit algebra(const KEY& k, const SCALAR& s = SCALAR(1))
        : VECT(k, s)
    {}

    /// Constructor from pointers to data range
    algebra(SCALAR const* begin, SCALAR const* end)
        : VECT(begin, end)
    {}

    /// Constructor from pointer to data range with offset
    algebra(DIMN offset, SCALAR const* begin, SCALAR const* end)
        : VECT(offset, begin, end)
    {}

    /// Constructor from pointer to data range with offset
    algebra(DIMN offset, SCALAR* begin, SCALAR* end)
        : VECT(offset, begin, end)
    {}

    template<typename InputIt, typename = typename std::iterator_traits<InputIt>::iterator_category>
    algebra(InputIt begin, InputIt end) : VECT(begin, end)
    {}

    template <typename B, typename M>
    explicit algebra(const algebra<B, Coeff, M, VectorType>& arg)
            : VECT(arg)
    {}

    algebra& operator=(const algebra&) = default;
    algebra& operator=(algebra&&) noexcept = default;

    const Multiplication& multiplication() const noexcept { return s_multiplication; }

public:

    //TODO: templatise these methods

    /// Adds to the instance a product of algebra instances.
    derived_type& add_mul(const algebra& a, const algebra& b)
    {
        auto& derived = static_cast<derived_type&>(*this);
        mtraits::multiply_and_add(s_multiplication, derived, a, b, scalar_passthrough());
        return derived;
    }

    /// Subtracts to the instance a product of algebra instances.
    derived_type& sub_mul(const algebra& a, const algebra& b)
    {
        auto& derived = static_cast<derived_type&>(*this);
        mtraits::multiply_and_add(s_multiplication, derived, a, b, scalar_minus());
        return derived;
    }

    /// Multiplies the instance by (algebra instance)*s.
    derived_type& mul_scal_prod(const algebra& rhs, const SCALAR& s)
    {
        auto& derived = static_cast<derived_type&>(*this);
        mtraits::multiply_inplace(s_multiplication, derived, rhs, scalar_post_mult(s));
        return derived;
    }

    /// Multiplies the instance by (algebra instance)*s up to maximum depth
    derived_type& mul_scal_prod(const algebra& rhs, const RATIONAL& s, const DEG depth)
    {
        auto& derived = static_cast<derived_type&>(*this);
        mtraits::mutiply_inplace(s_multiplication, derived, rhs, scalar_post_mult(s), depth);
        return derived;
    }

    /// Multiplies the instance by (algebra instance)/s.
    derived_type& mul_scal_div(const algebra& rhs, const RATIONAL& s)
    {
        auto& derived = static_cast<derived_type&>(*this);
        mtraits::multiply_inplace(s_multiplication, derived, rhs, rational_post_div(s));
        return derived;
    }

    /// Multiplies the instance by (algebra instance)/s up to maximum depth
    derived_type& mul_scal_div(const algebra& rhs, const RATIONAL& s, const DEG depth)
    {
        auto& derived = static_cast<derived_type&>(*this);
        mtraits::multiply_inplace(s_multiplication, derived, rhs, rational_post_div(s), depth);
        return derived;
    }

    /// Returns an instance of the commutator of two algebra instances.
    inline friend derived_type commutator(const algebra& a, const algebra& b)
    {// Returns a * b - b * a
        derived_type result;
        mtraits::multiply_and_add(s_multiplication, result, a, b, scalar_passthrough());
        mtraits::multiply_and_add(s_multiplication, result, b, a, scalar_minus());
        return result;
    }

    /// Returns a truncated version of the instance, by using basis::degree().
    inline derived_type truncate(const DEG min, const DEG max) const
    {
        derived_type result;
        const_iterator i;
        for (i = VECT::begin(); i != VECT::end(); ++i) {
            if ((VECT::basis.degree(i->key()) >= min) && (VECT::basis.degree(i->key()) <= max)) {
                result[i->key()] = i->value();
            }
        }
        return result;
    }

#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
private:
    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive& ar, unsigned int const /* version */)
    {
        ar& boost::serialization::base_object<VECT>(*this);
    }
#endif

public:
    template<typename Algebra1, typename OtherBasis, typename VType>
    friend
            typename std::enable_if<
                    std::is_base_of<algebra, Algebra1>::value && basis::related_to<OtherBasis, BASIS>::value,
                    Algebra1>::type
            operator*(const Algebra1& lhs, const vectors::dtl::vector<OtherBasis, Coeff, VType>& rhs)
    {
        Algebra1 result;
        mtraits::multiply_and_add(s_multiplication, result, lhs, rhs, scalar_passthrough());
        return result;
    }

    template<typename Algebra1, typename OtherBasis, typename VType>
    friend
            typename std::enable_if<
                    std::is_base_of<algebra, Algebra1>::value && basis::related_to<OtherBasis, BASIS>::value,
                    Algebra1>::type&
            operator*=(Algebra1& lhs, const vectors::dtl::vector<OtherBasis, Coeff, VType>& rhs)
    {
        mtraits::multiply_inplace(s_multiplication, lhs, rhs, scalar_passthrough());
        return lhs;
    }
};

template<typename B, typename C, typename M, template<typename, typename, typename...> class V, typename D, typename... Args>
const M algebra<B, C, M, V, D, Args...>::s_multiplication;


namespace dtl {

template <typename Algebra>
class is_algebra_impl
{
    template<typename B, typename C, typename M, template<typename, typename, typename...> class V, typename... Args>
    static std::true_type test(algebra<B, C, M, V, Args...>&);

    static std::false_type test(...);

public:
    static constexpr bool value = decltype(test(std::declval<Algebra&>()))::value;
};
} // namespace dtl


template <typename Algebra>
constexpr bool is_algebra() noexcept
{
    return dtl::is_algebra_impl<Algebra>::value;
}







}// namespace alg
// Include once wrapper
#endif// DJC_COROPA_LIBALGEBRA_ALGEBRAH_SEEN

// EOF.
