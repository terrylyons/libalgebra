//
// Created by sam on 12/02/2021.
//

#ifndef LIBALGEBRA_HYBRID_VECTOR_H
#define LIBALGEBRA_HYBRID_VECTOR_H

#include <algorithm>
#include <numeric>

#include "base_vector.h"
#include "dense_vector.h"
#include "detail/order_trait.h"
#include "sparse_vector.h"

#define DEFINE_FUSED_OP(NAME, ST, OP1, OP2)                            \
    hybrid_vector& NAME(const KEY& rhs, const ST s)                    \
    {                                                                  \
        DIMN idx;                                                      \
        if ((idx = DENSE::key_to_index(rhs)) < dense_dimension()) {    \
            DENSE::NAME(rhs, s);                                       \
        }                                                              \
        else {                                                         \
            SPARSE::NAME(rhs, s);                                      \
        }                                                              \
        maybe_resize();                                                \
        return *this;                                                  \
    }                                                                  \
                                                                       \
    hybrid_vector& NAME(const hybrid_vector& rhs, const ST s)          \
    {                                                                  \
        DIMN dim = std::max(dense_dimension(), rhs.dense_dimension()); \
        resize_dense_to_dimension(dim);                                \
        DENSE::NAME(rhs, s);                                           \
        SPARSE::NAME(rhs, s);                                          \
        maybe_resize();                                                \
        return *this;                                                  \
    }

namespace alg {
namespace vectors {

namespace tools {

struct size_control {

    template<typename Basis, typename Coeff, template<typename, typename> class Vector>
    static DIMN set_dense_dimension(vector<Basis, Coeff, Vector>& vect, DIMN dim)
    {
        Vector<Basis, Coeff>& v_vect = vect.base_vector();
        v_vect.resize_dense(dim);
        return v_vect.dense_dimension();
    }
};

}// namespace tools

namespace policy {

class basic_resize_policy
{
public:
    template<typename Vector, DEG MaxDegree>
    DIMN get_resize_size_impl(Vector const& vect, basis::with_degree<MaxDegree>)
    {
        DEG dense_deg = vect.dense_degree();
        if (dense_deg == MaxDegree) {
            return vect.max_dense_dimension();
        }

        std::vector<DIMN> degree_counts;
        degree_counts.resize(MaxDegree + 1);

        DEG d;
        typedef typename Vector::const_iterator citer;
        for (citer it(vect.sparse_begin()); it != vect.sparse_end(); ++it) {
            typename Vector::KEY const& key = it->key();
            d = vect.basis.degree(key);
            degree_counts[d] += 1;
        }

        DIMN degree_size;
        DEG resize_degree = dense_deg;
        for (d = dense_deg; d <= MaxDegree; ++d) {
            degree_size = vect.basis.start_of_degree(d + 1) - vect.basis.start_of_degree(d);
            if (degree_counts[d] >= (degree_size / 3)) {
                resize_degree = d;
            }
        }

        return vect.basis.start_of_degree(resize_degree + 1);
    }

    template<typename Vector>
    DIMN get_resize_size_impl(Vector const& vect, basis::without_degree)
    {
        DIMN dense_dim(vect.dense_dimension());
        DIMN sparse_dim(vect.sparse_size());

        auto info = basis::basis_traits<typename Vector::BASIS>::next_resize_dimension(vect.basis, dense_dim + 1, vect.dense_degree());
        DIMN next_dense_size(info.dimension);
        assert(next_dense_size <= vect.max_dense_dimension());
        assert(dense_dim <= vect.max_dense_dimension());

        if (sparse_dim > ((next_dense_size - dense_dim) / 4)) {
            return next_dense_size;
        }
        else {
            return dense_dim;
        }
    }

    template<typename Vector>
    DIMN get_resize_size(const Vector& vect)
    {
        return get_resize_size_impl(vect, vect.degree_tag);
    }
#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
private:
    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {}
#endif
};

}// namespace policy

namespace dtl {

#define LIBALGEBRA_HYBRID_REF_BINOP(OP)                              \
    template<typename T>                                             \
    hybrid_reference& operator OP(T arg)                             \
    {                                                                \
        if (m_is_dense) {                                            \
            p_vector->Hybrid::DENSE::operator[](m_key) OP scalar_type(arg);  \
        }                                                            \
        else {                                                       \
            p_vector->Hybrid::SPARSE::operator[](m_key) OP scalar_type(arg); \
        }                                                            \
        return *this;                                                \
    }

template<typename Hybrid>
class hybrid_reference
{
    using key_type = typename Hybrid::KEY;

    Hybrid* p_vector;
    const key_type& m_key;
    bool m_is_dense;

public:
    using scalar_type = typename Hybrid::SCALAR;

    explicit hybrid_reference(Hybrid* vector, const key_type& key, bool is_dense)
        : p_vector(vector), m_key(key), m_is_dense(is_dense)
    {}

    operator const scalar_type&() const noexcept
    {
        if (m_is_dense) {
            return static_cast<const scalar_type&>(p_vector->Hybrid::DENSE::operator[](m_key));
        }
        return static_cast<const scalar_type&>(p_vector->Hybrid::SPARSE::operator[](m_key));
    }

    LIBALGEBRA_HYBRID_REF_BINOP(=)
    LIBALGEBRA_HYBRID_REF_BINOP(+=)
    LIBALGEBRA_HYBRID_REF_BINOP(-=)
    LIBALGEBRA_HYBRID_REF_BINOP(*=)
    LIBALGEBRA_HYBRID_REF_BINOP(/=)
    LIBALGEBRA_HYBRID_REF_BINOP(&=)
    LIBALGEBRA_HYBRID_REF_BINOP(|=)
};

}// namespace dtl

/**
 * @brief Hybrid between a dense vector and a sparse vector.
 *
 * A hybrid vector is a combination of a dense vector for earlier keys in the basis and sparse for later elements.
 * The logic here is that the elements that appear earlier in the basis total order should be those which are
 * most used in various operations, while those that appear later are both less used and less likely to have
 * non-zero coefficients. When Basis has a degree, the dense part is resized according to the boundaries between
 * degrees, so that each degree level in its entirety is either held densely or sparsely.
 *
 * The dense part of the vector is dynamically resized, so that when a vector reaches a point at which a new index
 * or degree range is used frequently often, the dense part is resized to include this index/degree range. This
 * process is determined by the ResizePolicy template parameter.
 *
 * Operations are generally delegated separately to both parts of the vector, but some operations require
 * special logic to merge differently distributed data. In particular, the methods used to implement multiplication
 * for algebra types must.
 *
 * @tparam Basis Basis of the vector space
 * @tparam Coeffs Coefficient field
 * @tparam ResizePolicy Policy object used to determine when a vector should resize
 * @tparam SparseMap Type to use as the storage for the sparse vector
 */
template<typename Basis, typename Coeffs, typename ResizePolicy = policy::basic_resize_policy, typename SparseMap = std::unordered_map<typename Basis::KEY, typename Coeffs::S>>
class hybrid_vector : public dense_vector<Basis, Coeffs>, public sparse_vector<Basis, Coeffs, SparseMap>
{
    typedef dense_vector<Basis, Coeffs> DENSE;
    typedef sparse_vector<Basis, Coeffs, SparseMap> SPARSE;
    typedef ResizePolicy POLICY;

    friend struct tools::size_control;
    friend class dtl::data_access_base<hybrid_vector>;
    friend class dtl::hybrid_reference<hybrid_vector>;

private:
    // The resize manager is responsible for dictating when the
    // vector should resize it's dense part to contain a larger
    // proportion of the vector. This applies only to cases where
    // resizing is optional. This does not apply when resize is
    // mandatory.
    POLICY m_resize_policy;

public:
    // Type definitions
    typedef Basis BASIS;
    typedef typename BASIS::KEY KEY;

    typedef Coeffs COEFFS;
    typedef typename Coeffs::S SCALAR;
    typedef typename Coeffs::Q RATIONAL;

    typedef typename dtl::requires_order<Basis>::key_ordering key_ordering;

    using reference = dtl::hybrid_reference<hybrid_vector>;

public:
    // Static variables from base_Vec (via DENSE)
    using DENSE::basis;
    using DENSE::degree_tag;
    using DENSE::mone;
    using DENSE::one;
    using DENSE::zero;
    using typename SPARSE::coefficient_ring;

public:
    // Constructors

    /// Default constructor - creates the zero vector
    hybrid_vector()
        : DENSE(), SPARSE(), m_resize_policy()
    {}

    /// Create a new unidimensional vector with key and value
    explicit hybrid_vector(const KEY& key, const SCALAR& s = one)
        : DENSE(), SPARSE(key, s), m_resize_policy()
    {}

    /// Copy constructor
    hybrid_vector(const hybrid_vector& other)
        : DENSE(other), SPARSE(other), m_resize_policy()
    {}

    /// Move constructor
    hybrid_vector(hybrid_vector&& other) noexcept
        : DENSE(std::move(other.dense_part())), SPARSE(std::move(other.sparse_part())),
          m_resize_policy{std::move(other.m_resize_policy)}
    {}

    /**
     * @brief Construct from pointer to data
     *
     * Creates a new borrowed vector with dense elements from data
     *
     * @param begin Pointer to start of data
     * @param end Pointer to end of data
     */
    hybrid_vector(SCALAR const* begin, SCALAR const* end)
        : DENSE(begin, end), SPARSE()
    {}

    /**
     * @brief Construct a new denes vector with pointer to data and offset at beginning
     *
     * Creates a new dense vector where the dense elements are shifted by offset in the basis order
     *
     * @param offset size of offset before elements start
     * @param begin Pointer to start of data
     * @param end Pointer to end of data
     */
    hybrid_vector(DIMN offset, SCALAR const* begin, SCALAR const* end)
        : DENSE(offset, begin, end), SPARSE()
    {}

    /**
     * @brief Construct a new denes vector with pointer to data and offset at beginning
     *
     * Creates a new dense vector where the dense elements are shifted by offset in the basis order
     *
     * @param offset size of offset before elements start
     * @param begin Pointer to start of data
     * @param end Pointer to end of data
     */
    hybrid_vector(DIMN offset, SCALAR* begin, SCALAR* end)
        : DENSE(offset, begin, end), SPARSE()
    {}

    hybrid_vector& operator=(const hybrid_vector& other) = default;
    hybrid_vector& operator=(hybrid_vector&& other) noexcept = default;

private:
    /// Constructor from component vectors
    hybrid_vector(DENSE dense_vec, SPARSE sparse_vec)
        : DENSE(dense_vec), SPARSE(sparse_vec), m_resize_policy()
    {}

protected:
    // Resizing the dense part of the vector

    /// Reserve enough space in the dense part for given dimension
    void reserve_dense_to_dimension(const DIMN dim)
    {
        DENSE::resize_to_dimension(dim);
    }

    /// Reserve enough space in the dense part for elements of given degree
    void reserve_dense_to_degree(const DEG deg)
    {
        DENSE::resize_to_degree(deg);
    }

    /// Resize the dense part of the vector to contain the given dimension
    void resize_dense_to_dimension(const DIMN dim)
    {
        if (dim > dense_dimension()) {
            DENSE::resize_to_dimension(dim);
            incorporate_sparse();
        }
        else if (dim < dense_dimension()) {
            incorporate_dense(dim);
        }
    }

    /// Resize the dense part of a vector to contain the given degree
    void resize_dense_to_degree(const DEG deg)
    {
        resize_dense_to_dimension(DENSE::start_of_degree(deg));
    }

public:
    /// Trigger a resize and incorporate action if appropriate
    void maybe_resize()
    {
        DIMN resize_size = m_resize_policy.get_resize_size(*this);
        resize_dense_to_dimension(resize_size);
        incorporate_sparse();
    }

private:
    static bool sort_by_index(std::pair<DIMN, SCALAR> p1, std::pair<DIMN, SCALAR> p2)
    {
        return p1.first < p2.first;
    }

    /// Incorporate the sparse elements that should now be dense
    void incorporate_sparse()
    {
        DIMN dense_dim(dense_dimension());
        if (dense_dim == 0 || sparse_empty()) {
            return;
        }

        typename SPARSE::iterator it(SPARSE::begin()), end(SPARSE::end());
        std::vector<std::pair<DIMN, SCALAR>> buffer;
        buffer.reserve(sparse_size());

        DIMN idx;
        while (it != end) {
            if ((idx = key_to_index(it->key())) < dense_dim) {
                buffer.push_back(std::pair<DIMN, SCALAR>(idx, it->value()));
                // DENSE::value(idx) += it->value();
                SPARSE::erase(it++);
            }
            else {
                ++it;
            }
        }

        if (!utils::is_ordered<SparseMap>::value) {
            std::sort(buffer.begin(), buffer.end(), sort_by_index);
        }

        typename std::vector<std::pair<DIMN, SCALAR>>::const_iterator cit;
        for (cit = buffer.begin(); cit != buffer.end(); ++cit) {
            DENSE::value(cit->first) += cit->second;
        }
    }

    /// Incorporate the dense elements that will need to be sparse
    void incorporate_dense(const DIMN from_index)
    {
        if (from_index >= dense_dimension()) {
            return;
        }
        typedef std::pair<KEY, SCALAR> PAIR;
        std::vector<std::pair<KEY, SCALAR>> tmp;
        tmp.reserve(dense_dimension() - from_index);

        SCALAR val;
        for (DIMN i = from_index; i < dense_dimension(); ++i) {
            val = dense_value(i);
            if (val != zero) {
                tmp.push_back(PAIR(index_to_key(i), val));
            }
        }

        SPARSE::insert(tmp.begin(), tmp.end());
    }

public:
    // Vector information methods

    /// Get the dimension of the dense part of the vector
    DIMN dense_dimension() const
    {
        return DENSE::dimension();
    }

    /// Get the maximum possible dense dimension
    DIMN max_dense_dimension() const
    {
        return DENSE::max_dimension(degree_tag);
    }

    /// Get the maximum degree currently held by the dense part
    DEG dense_degree() const
    {
        return DENSE::degree();
    }

    /// Get the number of non-zero elements in the dense part
    DIMN dense_size() const
    {
        return DENSE::size();
    }

    /// Get the number of non-zero elements in the sparse part
    DIMN sparse_size() const
    {
        return SPARSE::size();
    }

    /// Get the number of non-zero elements in the vector
    DIMN size() const
    {
        return dense_size() + sparse_size();
    }


    DIMN dimension() const noexcept { return dense_dimension() + sparse_size(); }

    /// Test if the vector is empty (contains no non-zero elements)
    bool empty() const
    {
        return DENSE::empty() && SPARSE::empty();
    }

    /// Get the maximum degree of elements currently held by this vector
    DEG degree() const
    {
        return std::max(dense_degree(), SPARSE::degree());
    }

    /// Test if the maximum degree of elements in the vector is equal to given value
    bool degree_equals(const DEG degree) const
    {
        bool result(DENSE::degree_equals(degree));
        DEG d;
        for (typename SPARSE::const_iterator it(SPARSE::begin()); it != SPARSE::end(); ++it) {
            d = basis.degree(it->key());
            if (d > degree) {
                return false;
            }
            else if (!result && d == degree) {
                result = true;
            }
        }
        return result;
    }

    /// Test if the dense part is empty (contains no non-zero values)
    bool dense_empty() const
    {
        return DENSE::empty();
    }

    /// Test if the sparse part is empty (contains no non-zero values)
    bool sparse_empty() const
    {
        return SPARSE::empty();
    }

    // Sparse part and dense part access

    /// Get a reference to the dense part
    DENSE& dense_part()
    {
        return *this;
    }

    /// Get a const reference to the dense part
    const DENSE& dense_part() const
    {
        return *this;
    }

    /// get a reference to the sparse part
    SPARSE& sparse_part()
    {
        return *this;
    }

    /// Get a const reference to the sparse part
    const SPARSE& sparse_part() const
    {
        return *this;
    }

protected:
    // Information about keys from dense vector

    using DENSE::index_to_key;
    using DENSE::key_to_index;

public:
    class iterator_item
    {
        friend class iterators::vector_iterator<iterator_item>;

        friend class hybrid_vector;

    public:
        typedef KEY key_type;
        typedef SCALAR& value_type;

        iterator_item()
            : m_dense_iterator(), m_dense_end(), m_sparse_begin(), m_sparse_iterator()
        {}

        iterator_item(hybrid_vector& vect, typename SPARSE::iterator it)
            : m_dense_iterator(vect.dense_part().end()),
              m_dense_end(vect.dense_part().end()),
              m_sparse_begin(vect.sparse_part().begin()),
              m_sparse_iterator(it)
        {}

        iterator_item(hybrid_vector& vect, typename DENSE::iterator it)
            : m_dense_iterator(it),
              m_dense_end(vect.dense_part().end()),
              m_sparse_begin(vect.sparse_part().begin()),
              m_sparse_iterator(
                      vect.sparse_part().begin())
        {}

        key_type key() const
        {
            if (m_dense_iterator != m_dense_end) {
                return m_dense_iterator->key();
            }
            else {
                return m_sparse_iterator->key();
            }
        }

        value_type value() const
        {
            if (m_dense_iterator != m_dense_end) {
                return m_dense_iterator->value();
            }
            else {
                return m_sparse_iterator->value();
            }
        }

        DIMN index() const
        {
            if (is_dense()) {
                return m_dense_iterator->index();
            }
            else {
                return m_dense_iterator->index();
            }
        }

    private:
        typename DENSE::iterator m_dense_iterator;
        typename DENSE::iterator m_dense_end;
        typename SPARSE::iterator m_sparse_begin;
        typename SPARSE::iterator m_sparse_iterator;

    private:
        bool compare_iterators(const iterator_item& other) const
        {
            assert(m_sparse_begin == other.m_sparse_begin);
            assert(m_dense_end == other.m_dense_end);

            bool lhs_dense(is_dense()), rhs_dense(other.is_dense());
            if (lhs_dense && rhs_dense) {
                return m_dense_iterator == other.m_dense_iterator;
            }
            else if (lhs_dense || rhs_dense) {
                return false;
            }
            else {
                return m_sparse_iterator == other.m_sparse_iterator;
            }
        }

        void advance()
        {
            if (m_dense_iterator != m_dense_end) {
                ++m_dense_iterator;
            }
            else {
                ++m_sparse_iterator;
            }
        }

        bool is_dense() const
        {
            return m_dense_iterator != m_dense_end;
        }
    };

    class const_iterator_item
    {
        friend class iterators::vector_iterator<const_iterator_item>;

        friend class hybrid_vector;

    public:
        typedef KEY key_type;
        typedef const SCALAR& value_type;

        const_iterator_item()
            : m_dense_iterator(), m_dense_end(), m_sparse_begin(), m_sparse_iterator()
        {}

        const_iterator_item(const hybrid_vector& vect, typename SPARSE::const_iterator it)
            : m_dense_iterator(
                    vect.dense_part().end()),
              m_dense_end(vect.dense_part().end()), m_sparse_begin(
                                                            vect.sparse_part().begin()),
              m_sparse_iterator(it)
        {}

        const_iterator_item(const hybrid_vector& vect, typename DENSE::const_iterator it)
            : m_dense_iterator(it),
              m_dense_end(
                      vect.dense_part()
                              .end()),
              m_sparse_begin(
                      vect.sparse_part()
                              .begin()),
              m_sparse_iterator(
                      vect.sparse_part()
                              .begin())
        {
            assert(vect.sparse_empty() || m_sparse_begin != vect.sparse_part().end());
        }

        key_type key() const
        {
            if (is_dense()) {
                return m_dense_iterator->key();
            }
            else {
                return m_sparse_iterator->key();
            }
        }

        value_type value() const
        {
            if (is_dense()) {
                return m_dense_iterator->value();
            }
            else {
                return m_sparse_iterator->value();
            }
        }

        DIMN index() const
        {
            if (is_dense()) {
                return m_dense_iterator->index();
            }
            else {
                return m_dense_iterator->index();
            }
        }

    private:
        typename DENSE::const_iterator m_dense_iterator;
        typename DENSE::const_iterator m_dense_end;
        typename SPARSE::const_iterator m_sparse_begin;
        typename SPARSE::const_iterator m_sparse_iterator;

    private:
        bool compare_iterators(const const_iterator_item& other) const
        {
            assert(m_sparse_begin == other.m_sparse_begin);
            assert(m_dense_end == other.m_dense_end);

            bool lhs_dense(is_dense()), rhs_dense(other.is_dense());
            if (lhs_dense && rhs_dense) {
                return m_dense_iterator == other.m_dense_iterator;
            }
            else if (lhs_dense || rhs_dense) {
                return false;
            }
            else {
                return m_sparse_iterator == other.m_sparse_iterator;
            }
        }

        void advance()
        {
            if (is_dense()) {
                ++m_dense_iterator;
            }
            else {
                ++m_sparse_iterator;
            }
        }

        bool is_dense() const
        {
            return m_dense_iterator != m_dense_end;
        }
    };

    typedef iterators::vector_iterator<iterator_item> iterator;
    typedef iterators::vector_iterator<const_iterator_item> const_iterator;

    // Iterator methods

    /// Iterator to beginning of vector
    iterator begin()
    {
        if (dense_dimension() == 0) {
            return iterator(*this, SPARSE::begin());
        }
        else {
            return iterator(*this, DENSE::begin());
        }
    }

    /// Iterator to end of vector
    iterator end()
    {
        return iterator(*this, SPARSE::end());
    }

    /// Const iterator to beginning of vector
    const_iterator begin() const
    {
        if (dense_dimension() == 0) {
            return const_iterator(*this, SPARSE::begin());
        }
        else {
            return const_iterator(*this, DENSE::begin());
        }
    }

    /// Const iterator to end of vector
    const_iterator end() const
    {
        return const_iterator(*this, SPARSE::end());
    }

    /// Const iterator to beginning of vector
    const_iterator cbegin() const
    {
        return begin();
    }

    /// Const iterator to end of vector
    const_iterator cend() const
    {
        return end();
    }

    /// Const iterator to beginning of densely stored elements
    const_iterator sparse_begin() const
    {
        return const_iterator(*this, SPARSE::begin());
    }
    /// Const iterator to end of sparsely stored elements
    const_iterator sparse_end() const
    {
        return const_iterator(*this, SPARSE::end());
    }

    /// Insert element with iterator with hint iterator
    iterator insert(iterator it, const std::pair<const KEY, SCALAR>& val)
    {
        std::pair<iterator, bool> rv = insert(val);
        return rv.first;
    }

    /// Insert an element into the vector
    std::pair<iterator, bool> insert(const KEY key, SCALAR val)
    {
        std::pair<iterator, bool> rv(end(), false);
        DIMN idx;
        if ((idx = key_to_index(key)) < dense_dimension()) {
            if (DENSE::value(idx) == zero) {
                DENSE::value(idx) = val;
                rv.second = true;
            }
            rv.first = iterator(*this, DENSE::find_index(idx));
        }
        else {
            std::pair<const KEY, SCALAR> tmp_pair(key, val);
            std::pair<typename SPARSE::iterator, bool> tmp = SPARSE::insert(tmp_pair);
            rv.first = iterator(*this, tmp.first);
            rv.second = tmp.second;
        }
        return rv;
    }

    /// Insert an element into the vector from key-value pair
    std::pair<iterator, bool> insert(const std::pair<const KEY, SCALAR>& p)
    {
        DIMN idx;
        if ((idx = key_to_index(p.first)) < dense_dimension()) {
            if (DENSE::value(idx) == zero && p.second != zero) {
                DENSE::value(idx) = p.second;
                return std::pair<iterator, bool>(iterator(*this, DENSE::find_index(idx)), true);
            }
        }
        else {
            //typename SPARSE::iterator it = SPARSE::find(p.first);
            if (p.second != zero) {
                std::pair<typename SPARSE::iterator, bool> ins = SPARSE::insert(p);
                return std::pair<iterator, bool>(iterator(*this, ins.first), ins.second);
            }
        }
        return std::pair<iterator, bool>(end(), false);
    }

    /// Insert elements from iterator
    template<typename InputIterator>
    void insert(InputIterator begin, InputIterator end)
    {
        for (InputIterator it(begin); it != end; ++it) {
            insert(it->first, it->second);
        }
    }

    /// Find the const iterator of a particular key in the vector
    const_iterator find(const KEY& key) const
    {
        DIMN idx;
        if ((idx = key_to_index(key)) < dense_dimension()) {
            return const_iterator(*this, DENSE::find_index(idx));
        }
        else {
            return const_iterator(*this, SPARSE::find(key));
        }
    }

    /// Find the iterator of a particular key in the vector
    iterator find(const KEY& key)
    {
        DIMN idx;
        if ((idx = key_to_index(key)) < dense_dimension()) {
            return iterator(*this, DENSE::find_index(idx));
        }
        else {
            return iterator(*this, SPARSE::find(key));
        }
    }

    /// Erase (set to zero) the coefficient of given key
    void erase(const KEY& key)
    {
        DIMN idx;
        if ((idx = key_to_index(key)) < dense_dimension()) {
            DENSE::value(idx) = zero;
        }
        else {
            SPARSE::erase(key);
        }
    }

    /// Erase (set to zero) the value associated with an iterator
    void erase(iterator& it)
    {
        if (it->is_dense()) {
            it->value() = zero;
        }
        else {
            SPARSE::erase(it->m_sparse_iterator);
        }
    }

public:
    /// Clear the vector (set all elements to zero)
    void clear()
    {
        DENSE::clear();
        SPARSE::clear();
    }

public:
    /// Swap the data held in this vector with another
    void swap(hybrid_vector& other)
    {
        DENSE::swap(other);
        SPARSE::swap(other);
        maybe_resize();
    }

public:
    // Element access

    /// Get const reference to coefficient of Key
    const SCALAR& operator[](const KEY& key) const
    {
        DIMN idx;
        if ((idx = key_to_index(key)) < DENSE::dimension()) {
            return DENSE::value(idx);
        }
        else {
            return SPARSE::operator[](key);
        }
    }

    /// Get reference to coefficient of key
    reference operator[](const KEY& key)
    {
        DIMN idx = key_to_index(key);
        return reference(this, key, idx < dense_dimension());
    }

    /// Get a reference to the coefficient of the key using index
    SCALAR& value(const DIMN idx)
    {
        if (idx < DENSE::dimension()) {
            return DENSE::value(idx);
        }
        else {
            return SPARSE::operator[](index_to_key(idx));
        }
    }

    /// Get a const reference to the coefficient of the key using index
    const SCALAR& value(const DIMN idx) const
    {
        if (idx < DENSE::dimension()) {
            return DENSE::value(idx);
        }
        else {
            return SPARSE::operator[](index_to_key(idx));
        }
    }

    //    /// Get reference to an element from the dense part at index
    //    SCALAR& dense_value(const DIMN idx)
    //    {
    //        return DENSE::value(idx);
    //    }

    /// Get a const reference to an element from the dense part at index
    const SCALAR& dense_value(const DIMN idx) const
    {
        return DENSE::value(idx);
    }

public:
    // Comparison operators

    /// Equality operator
    bool operator==(const hybrid_vector& rhs) const
    {
        if (dense_dimension() == 0 && rhs.dense_dimension() == 0) {
            return sparse_part() == rhs.sparse_part();
        }

        if (sparse_empty() && rhs.sparse_empty()) {
            return dense_part() == rhs.dense_part();
        }

        std::pair<DIMN, bool> dense_part_eq = DENSE::equal_to_min(rhs);
        if (dense_part_eq.first != 0 && !dense_part_eq.second) {
            return false;
        }

        typename SPARSE::const_iterator cit, cend, oit;

        if (dense_part_eq.first == dense_dimension() && dense_part_eq.first == rhs.dense_dimension()) {

            return sparse_part() == rhs.sparse_part();
        }
        else if (dense_part_eq.first == dense_dimension()) {
            // Case rhs.dense_dimension() > dense_dimension()

            cend = SPARSE::end();
            DIMN sz(0), sparse_sz(sparse_size());

            for (DIMN i = dense_part_eq.first; i < rhs.dense_dimension(); ++i) {
                if (rhs.dense_value(i) == zero) {
                    continue;
                }

                cit = SPARSE::find(index_to_key(i));
                if (cit != cend) {
                    if (cit->value() != rhs.dense_part().value(i)) {
                        return false;
                    }
                    sz += 1;
                }
                else {
                    return false;
                }
            }

            for (oit = rhs.sparse_part().begin(); oit != rhs.sparse_part().end(); ++oit) {
                cit = SPARSE::find(oit->key());
                if (cit != cend) {
                    if (cit->value() != oit->value()) {
                        return false;
                    }
                    sz += 1;
                }
                else {
                    return false;
                }
            }

            return sz == sparse_sz;
        }
        else if (dense_part_eq.first == rhs.dense_dimension()) {
            cend = rhs.sparse_part().end();

            cend = SPARSE::end();
            DIMN sz(0), sparse_sz(rhs.sparse_size());

            for (DIMN i = dense_part_eq.first; i < dense_dimension(); ++i) {
                if (dense_value(i) == zero) {
                    continue;
                }

                cit = rhs.sparse_part().find(index_to_key(i));
                if (cit != cend) {
                    if (cit->value() != rhs.dense_part().value(i)) {
                        return false;
                    }
                    sz += 1;
                }
                else {
                    return false;
                }
            }

            for (oit = sparse_part().begin(); oit != sparse_part().end(); ++oit) {
                cit = rhs.sparse_part().find(oit->key());
                if (cit != cend) {
                    if (cit->value() != oit->value()) {
                        return false;
                    }
                    sz += 1;
                }
                else {
                    return false;
                }
            }

            return sz == sparse_sz;
        }

        // Should be unreachable
        assert(false);
        return false;
    }

    /// Non-equality operator
    bool operator!=(const hybrid_vector& other) const
    {
        return !operator==(other);
    }

public:
    template<typename F>
    static void apply_unary_operation(
            hybrid_vector& result,
            const hybrid_vector& arg,
            F&& func)
    {
        DENSE::apply_unary_operation(result, arg, func);
        SPARSE::apply_unary_operation(result, arg, func);
    }

    template<typename F>
    static void apply_inplace_unary_op(hybrid_vector& arg, F&& func)
    {
        DENSE::apply_inplace_unary_op(arg, func);
        SPARSE::apply_inplace_unary_op(arg, func);
    }

    template<typename F>
    static void apply_flat_binary_operation(
            hybrid_vector& result,
            const hybrid_vector& lhs,
            const hybrid_vector& rhs,
            F&& func)
    {
        DENSE::apply_flat_binary_operation(result, lhs, rhs, func);
        SPARSE::apply_flat_binary_operation(result, lhs, rhs, func);
        result.incorporate_sparse();
    }

    template<typename F>
    static void apply_inplace_flat_binary_op(
            hybrid_vector& lhs,
            const hybrid_vector& rhs,
            F&& func)
    {
        DENSE::apply_inplace_flat_binary_op(lhs, rhs, func);
        SPARSE::apply_inplace_flat_binary_op(lhs, rhs, func);
        lhs.incorporate_sparse();
    }

    // Arithmetic operators
    // The first few are operators that are simply applied to each
    // part individually.

    /// Unary minus operation
    hybrid_vector operator-(void) const
    {
        return hybrid_vector(DENSE::operator-(), SPARSE::operator-());
    }

    /// Scalar multiplication
    hybrid_vector& operator*=(const SCALAR s)
    {
        DENSE::operator*=(s);
        SPARSE::operator*=(s);
        return *this;
    }

    /// Rational division
    hybrid_vector& operator/=(const RATIONAL s)
    {
        DENSE::operator/=(s);
        SPARSE::operator/=(s);
        return *this;
    }

    // The next few operators require some matching
    /// Inplace addition
    hybrid_vector& operator+=(const hybrid_vector& other)
    {
        DIMN dim = std::max(dense_dimension(), other.dense_dimension());
        resize_dense_to_dimension(dim);

        DENSE::operator+=(other);
        SPARSE::operator+=(other);
        maybe_resize();
        return *this;
    }

    /// Inplace subtraction
    hybrid_vector& operator-=(const hybrid_vector& other)
    {
        DIMN dim = std::max(dense_dimension(), other.dense_dimension());
        resize_dense_to_dimension(dim);

        DENSE::operator-=(other);
        SPARSE::operator-=(other);
        maybe_resize();
        return *this;
    }

    /// Inplace coordinatewise minimum
    hybrid_vector& operator&=(const hybrid_vector& rhs)
    {
        DIMN dim = std::max(dense_dimension(), rhs.dense_dimension());
        resize_dense_to_dimension(dim);

        DENSE::operator&=(rhs);
        SPARSE::operator&=(rhs);

        maybe_resize();
        return *this;
    }

    /// Inplace coordinatewise minimum
    hybrid_vector& operator|=(const hybrid_vector& rhs)
    {
        DIMN dim = std::max(dense_dimension(), rhs.dense_dimension());
        resize_dense_to_dimension(dim);

        DENSE::operator|=(rhs);
        SPARSE::operator|=(rhs);

        maybe_resize();
        return *this;
    }

    /// Scalar multiplication
    hybrid_vector operator*(const SCALAR& rhs) const
    {
        hybrid_vector result(*this);
        return result *= rhs;
    };
    /// Rational division
    hybrid_vector operator/(const RATIONAL& rhs) const
    {
        hybrid_vector result(*this);
        return result /= rhs;
    };
    /// Addition
    hybrid_vector operator+(const hybrid_vector& rhs) const
    {
        hybrid_vector result(*this);
        return result += rhs;
    };
    /// Subtraction
    hybrid_vector operator-(const hybrid_vector& rhs) const
    {
        hybrid_vector result(*this);
        return result -= rhs;
    };
    /// Coordinatewise minimum
    hybrid_vector operator&(const hybrid_vector& rhs) const
    {
        hybrid_vector result(*this);
        return result &= rhs;
    };
    /// Coordinatewise maximum
    hybrid_vector operator|(const hybrid_vector& rhs) const
    {
        hybrid_vector result(*this);
        return result |= rhs;
    };

public:
    // Fused operations
    /// Fused inplace addition with scalar multiply
    DEFINE_FUSED_OP(add_scal_prod, SCALAR, +=, *);
    /// Fused inplace subtraction with scalar multiply
    DEFINE_FUSED_OP(sub_scal_prod, SCALAR, -=, *);
    /// Fused inplace addition with rational divide
    DEFINE_FUSED_OP(add_scal_div, RATIONAL, +=, /);
    /// Fused inplace subtraction with rational divide
    DEFINE_FUSED_OP(sub_scal_div, RATIONAL, -=, /);

public:
    // Norms

    /// Compute the L1 norm of the vector
    SCALAR NormL1() const
    {
        return DENSE::NormL1() + SPARSE::NormL1();
    }

    /// Compute the L1 norm of elements up to given degree
    SCALAR NormL1(const DEG deg) const
    {
        return DENSE::NormL1(deg) + SPARSE::NormL1(deg);
    }

    /// Compute the L-infinity norm of the vector
    SCALAR NormLInf() const
    {
        SCALAR dli = DENSE::NormLInf();
        SCALAR sli = SPARSE::NormLInf();
        return std::max(dli, sli);
    }

    /// Compute the L-infinity norm of elements up to maximum degree
    SCALAR NormLInf(const DEG deg) const
    {
        return std::max(DENSE::NormLInf(deg), SPARSE::NormLInf(deg));
    }

protected:
    // Display

    /// Print members to output stream
    void print_members(std::ostream& os) const
    {
        DENSE::print_members(os);
        SPARSE::print_members(os);
    }

public:
    /// Output stream operator
    inline friend std::ostream& operator<<(std::ostream& os, const hybrid_vector& rhs)
    {
        os << '{';
        rhs.print_members(os);
        os << ' ' << '}';
        return os;
    }

public:
    /// Comparison operator for key-value pairs
    static bool comp(std::pair<KEY, SCALAR> lhs, std::pair<KEY, SCALAR> rhs)
    {
        key_ordering ord;
        return ord(lhs.first, rhs.first);
    }

#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
private:
    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& boost::serialization::base_object<DENSE>(*this);
        ar& boost::serialization::base_object<SPARSE>(*this);
        ar& m_resize_policy;
    }
#endif
};

#undef DEFINE_FUSED_OP

namespace dtl {

template<typename Basis, typename Coeffs, typename Policy, typename Map>
struct data_access_base<hybrid_vector<Basis, Coeffs, Policy, Map>> {

    using vector_type = hybrid_vector<Basis, Coeffs, Policy, Map>;
    using tag = access_type_sparse;

    static typename vector_type::const_iterator range_begin(const vector_type& vect)
    {
        return vect.begin();
    }

    static typename vector_type::const_iterator range_end(const vector_type& vect)
    {
        return vect.end();
    }

    static typename vector_type::iterator range_begin(vector_type& vect)
    {
        return vect.begin();
    }

    static typename vector_type::iterator range_end(vector_type& vect)
    {
        return vect.end();
    }
};

}// namespace dtl

}// namespace vectors
}// namespace alg

#endif// LIBALGEBRA_HYBRID_VECTOR_H
