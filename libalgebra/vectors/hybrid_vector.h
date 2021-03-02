//
// Created by sam on 12/02/2021.
//

#ifndef LIBALGEBRA_HYBRID_VECTOR_H
#define LIBALGEBRA_HYBRID_VECTOR_H

#include "libalgebra/vectors/base_vector.h"
#include "libalgebra/vectors/sparse_vector.h"
#include "libalgebra/vectors/dense_vector.h"


#define DEFINE_FUSED_OP(NAME, ST, OP1, OP2)                                 \
    hybrid_vector& NAME(const KEY& rhs, const ST s)                         \
    {                                                                       \
        DIMN idx;                                                           \
        if ((idx = DENSE::key_to_index(rhs)) < DENSE::dimension()) {        \
            DENSE::value(idx) OP1 (one OP2 s);                              \
        } else {                                                            \
            SPARSE::operator[](rhs) OP1 (one OP2 s);                        \
        }                                                                   \
        return *this;                                                       \
    }                                                                       \
                                                                            \
    hybrid_vector& NAME(const hybrid_vector& rhs, const ST s)               \
    {                                                                       \
        DIMN dim = std::max(DENSE::dimension(), rhs.dimension());           \
        resize_dense_to_dimension(dim);                                     \
        DENSE:: NAME (rhs, s);                                              \
        SPARSE:: NAME (rhs, s);                                             \
        incorporate_sparse();                                               \
        return *this;                                                       \
    }

namespace alg {
namespace vectors {


namespace tools {

struct size_control
{

    template <typename Basis, typename Coeff, typename Vector>
    static DIMN set_dense_dimension(vector<Basis, Coeff, Vector>& vect, DIMN dim)
    {
        Vector& v_vect = dtl::vector_base_access::convert(vect);
        v_vect.resize_dense(dim);
        return v_vect.dense_dimension();
    }


};

}

class basic_resize_manager
{
public:

};


template<
        typename Basis,
        typename Coeffs,
        typename ResizeManager=basic_resize_manager,
        typename DenseStorage=std::vector<typename Coeffs::S>,
        typename SparseMap=LIBALGEBRA_DEFAULT_MAP_TYPE >
class hybrid_vector : public dense_vector<Basis, Coeffs, DenseStorage>,
                      public sparse_vector<Basis, Coeffs, SparseMap>
{
    typedef dense_vector <Basis, Coeffs, DenseStorage> DENSE;
    typedef sparse_vector <Basis, Coeffs, SparseMap> SPARSE;
    typedef ResizeManager MANAGER;

    friend struct tools::size_control;

private:

    // The resize manager is responsible for dictating when the
    // vector should resize it's dense part to contain a larger
    // proportion of the vector. This applies only to cases where
    // resizing is optional. This does not apply when resize is
    // mandatory.
    MANAGER m_resize_manager;

public:

    // Type definitions
    typedef Basis BASIS;
    typedef typename BASIS::KEY KEY;

    typedef typename Coeffs::S SCALAR;
    typedef typename Coeffs::Q RATIONAL;

    typedef typename dtl::requires_order<Basis>::key_ordering key_ordering;


public:

    // Static variables from base_Vec (via DENSE)
    using DENSE::basis;
    using DENSE::mone;
    using DENSE::zero;
    using DENSE::one;
    using DENSE::degree_tag;

public:

    // Constructors

    hybrid_vector(void) : DENSE(), SPARSE()
    {}

    explicit hybrid_vector(const KEY &key, const SCALAR& s = one)
            : DENSE(), SPARSE(key, s)
    {
    }

    hybrid_vector(const hybrid_vector &other)
            : DENSE(other), SPARSE(other)
    {}

private:

    hybrid_vector(DENSE dense_vec, SPARSE sparse_vec)
            : DENSE(dense_vec), SPARSE(sparse_vec)
    {
    }

protected:

    // Resizing the dense part of the vector

    void reserve_dense_to_dimension(const DIMN dim)
    {
        DENSE::resize_to_dimension(dim);
    }

    void reserve_dense_to_degree(const DEG deg)
    {
        DENSE::resize_to_degree(deg);
    }

    void resize_dense_to_dimension(const DIMN dim)
    {
        if (dim > dense_dimension()) {
            DENSE::resize_to_dimension(dim);
            incorporate_sparse();
        } else if (dim < dense_dimension()) {
            incorporate_dense(dim);
            DENSE::resize_to_dimension(dim);
        }
    }

    void resize_dense_to_degree(const DEG deg)
    {
        resize_dense_to_dimension(DENSE::start_of_degree(deg));
    }


    void maybe_resize()
    {
        DIMN sparse_sz(sparse_size());
        while (sparse_sz > dense_dimension()) {
            resize_dense_to_dimension(2*dense_dimension() + 1);
        }
        incorporate_sparse();
    }

private:

    /// Incorporate the sparse elements that should now be dense
    void incorporate_sparse()
    {
        if (dense_dimension() == 0 || sparse_empty()) {
            return;
        }

        typename SPARSE::iterator it(SPARSE::begin()), end(SPARSE::end());


        DIMN idx;

        while (it != end) {
            if ((idx = key_to_index(it->key())) < dense_dimension()) {
                DENSE::value(idx) += it->value();
                SPARSE::erase(it++);
            } else {
                ++it;
            }
        }

    }

    /// Incorporate the dense elements that will need to be sparse
    void incorporate_dense(const DIMN from_index)
    {
        if (from_index >= dense_dimension()) {
            return;
        }
        typedef std::pair<KEY, SCALAR> PAIR;
        std::vector<std::pair<KEY, SCALAR> > tmp;
        tmp.reserve(dense_dimension() - from_index);

        for (DIMN i = from_index; i < dense_dimension(); ++i) {
            tmp.push_back(PAIR(index_to_key(i), DENSE::value(i)));
        }

        DENSE::resize_to_dimension(from_index);
        SPARSE::insert(tmp.begin(), tmp.end());

    }

public:

    // Vector information methods

    DIMN dense_dimension() const
    { return DENSE::dimension(); }

    DIMN max_dense_dimension() const
    {
        return DENSE::max_dimension(degree_tag);
    }

    DEG dense_degree() const
    { return DENSE::degree(); }

    DIMN dense_size() const
    { return DENSE::size(); }

    DIMN sparse_size() const
    { return SPARSE::size(); }

    DIMN size() const
    {
        return dense_size() + sparse_size();
    }

    bool empty() const
    {
        return DENSE::empty() && SPARSE::empty();
    }

    DEG degree() const
    {
        return std::max(dense_degree(), SPARSE::degree());
    }

    bool degree_equals(const DEG degree) const
    {
        bool result (DENSE::degree_equals(degree));
        DEG d;
        for (typename SPARSE::const_iterator it(SPARSE::begin()); it != SPARSE::end(); ++it) {
            d = basis.degree(it->key());
            if (d > degree) {
                return false;
            } else if (!result && d == degree) {
                result = true;
            }
        }
        return result;
    }


    bool dense_empty() const
    { return DENSE::empty(); }

    bool sparse_empty() const
    { return SPARSE::empty(); }

protected:

    // Sparse part and dense part access

    DENSE &dense_part()
    { return *this; }

    const DENSE &dense_part() const
    { return *this; }

    SPARSE &sparse_part()
    { return *this; }

    const SPARSE &sparse_part() const
    { return *this; }

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
        typedef SCALAR &value_type;

        iterator_item()
                : m_dense_iterator(),
                  m_dense_end(),
                  m_sparse_begin(),
                  m_sparse_iterator()
        {}

        iterator_item(const iterator_item &other)
                : m_dense_iterator(other.m_dense_iterator),
                  m_dense_end(other.m_dense_end),
                  m_sparse_begin(other.m_sparse_begin),
                  m_sparse_iterator(other.m_sparse_iterator)
        {}

        iterator_item(hybrid_vector &vect, typename SPARSE::iterator it)
                : m_dense_iterator(vect.dense_part().end()),
                  m_dense_end(vect.dense_part().end()),
                  m_sparse_begin(vect.sparse_part().begin()),
                  m_sparse_iterator(it)
        {}

        iterator_item(hybrid_vector &vect, typename DENSE::iterator it)
                : m_dense_iterator(it),
                  m_dense_end(vect.dense_part().end()),
                  m_sparse_begin(vect.sparse_part().begin()),
                  m_sparse_iterator(m_sparse_begin)
        {}

        key_type key()
        {
            if (m_dense_iterator != m_dense_end) {
                return m_dense_iterator->key();
            } else {
                return m_sparse_iterator->key();
            }
        }

        value_type value()
        {
            if (m_dense_iterator != m_dense_end) {
                return m_dense_iterator->value();
            } else {
                return m_sparse_iterator->value();
            }
        }


    private:

        typename DENSE::iterator m_dense_iterator;
        typename DENSE::iterator m_dense_end;
        typename SPARSE::iterator m_sparse_iterator;
        typename SPARSE::iterator m_sparse_begin;

    private:

        bool compare_iterators(const iterator_item &other) const
        {
            assert(m_sparse_begin == other.m_sparse_begin);
            assert(m_dense_end == other.m_dense_end);
            if (is_dense() && other.is_dense()) {
                return m_dense_iterator == other.m_dense_iterator;
            } else {
                return m_sparse_iterator == other.m_sparse_iterator;
            }
        }

        void advance()
        {
            if (m_dense_iterator != m_dense_end) {
                ++m_dense_iterator;
            } else {
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
        typedef const SCALAR &value_type;

        const_iterator_item()
                : m_dense_iterator(),
                  m_dense_end(),
                  m_sparse_begin(),
                  m_sparse_iterator()
        {}

        const_iterator_item(const const_iterator_item &other)
                : m_dense_iterator(other.m_dense_iterator),
                  m_dense_end(other.m_dense_end),
                  m_sparse_begin(other.m_sparse_begin),
                  m_sparse_iterator(other.m_sparse_iterator)
        {}

        const_iterator_item(const hybrid_vector &vect, typename SPARSE::const_iterator it)
                : m_dense_iterator(vect.dense_part().end()),
                  m_dense_end(vect.dense_part().end()),
                  m_sparse_begin(vect.sparse_part().begin()),
                  m_sparse_iterator(it)
        {}

        const_iterator_item(const hybrid_vector &vect, typename DENSE::const_iterator it)
                : m_dense_iterator(it),
                  m_dense_end(vect.dense_part().end()),
                  m_sparse_begin(vect.sparse_part().begin()),
                  m_sparse_iterator(m_sparse_begin)
        {
            assert(vect.sparse_empty() || m_sparse_begin != vect.sparse_part().end());
        }

        key_type key()
        {
            if (m_dense_iterator != m_dense_end) {
                return m_dense_iterator->key();
            } else {
                return m_sparse_iterator->key();
            }
        }

        value_type value()
        {
            if (m_dense_iterator != m_dense_end) {
                return m_dense_iterator->value();
            } else {
                return m_sparse_iterator->value();
            }
        }


    private:

        typename DENSE::const_iterator m_dense_iterator;
        typename DENSE::const_iterator m_dense_end;
        typename SPARSE::const_iterator m_sparse_iterator;
        typename SPARSE::const_iterator m_sparse_begin;

    private:

        bool compare_iterators(const const_iterator_item &other) const
        {
            assert(m_sparse_begin == other.m_sparse_begin);
            assert(m_dense_end == other.m_dense_end);
            if (is_dense() && other.is_dense()) {
                return m_dense_iterator == other.m_dense_iterator;
            } else {
                return m_sparse_iterator == other.m_sparse_iterator;
            }
        }

        void advance()
        {
            if (m_dense_iterator != m_dense_end) {
                ++m_dense_iterator;
            } else {
                ++m_sparse_iterator;
            }
        }

        bool is_dense() const
        {
            return m_dense_iterator != m_dense_end;
        }

    };

    typedef iterators::vector_iterator <iterator_item> iterator;
    typedef iterators::vector_iterator <const_iterator_item> const_iterator;

    // Iterator methods

    iterator begin()
    {
        if (dense_empty()) {
            return iterator(*this, SPARSE::begin());
        } else {
            return iterator(*this, DENSE::begin());
        }
    }

    iterator end()
    {
        return iterator(*this, SPARSE::end());
    }

    const_iterator begin() const
    {
        if (dense_empty()) {
            return const_iterator(*this, SPARSE::begin());
        } else {
            return const_iterator(*this, DENSE::begin());
        }
    }

    const_iterator end() const
    {
        return const_iterator(*this, SPARSE::end());
    }

    const_iterator cbegin() const
    {
        return begin();
    }

    const_iterator cend() const
    {
        return end();
    }

    iterator insert(iterator it, const std::pair<const KEY, SCALAR> &val)
    {
        std::pair<iterator, bool> rv = insert(val);
        return rv.first;
    }

    std::pair<iterator, bool> insert(const KEY key, SCALAR val)
    {
        std::pair<iterator, bool> rv(end(), false);
        DIMN idx;
        if ((idx = key_to_index(key)) < dense_dimension()) {
            if (DENSE::value(idx) == zero) {
                DENSE::value(idx) = val;
                rv.second = true;
            }
            rv.first = iterator(*this, DENSE::find(idx));
        } else {
            std::pair<const KEY, SCALAR> tmp_pair(key, val);
            std::pair<typename SPARSE::iterator, bool> tmp =
                    SPARSE::insert(tmp_pair);
            rv.first = iterator(*this, tmp.first);
            rv.second = tmp.second;
        }
        return rv;
    }

    std::pair<iterator, bool> insert(const std::pair<const KEY, SCALAR> &p)
    {
        DIMN idx;
        if ((idx = key_to_index(p.first)) < dense_dimension()) {
            if (DENSE::value(idx) == zero && p.second != zero) {
                DENSE::value(idx) = p.second;
                return std::pair<iterator, bool>(
                        iterator(*this, DENSE::find(idx)),
                        true);
            }
        } else {
            typename SPARSE::iterator it = SPARSE::find(p.first);
            if (p.second != zero) {
                std::pair<typename SPARSE::iterator, bool>
                        ins = SPARSE::insert(p);
                return std::pair<iterator, bool>(
                        iterator(*this, ins.first),
                        ins.second);
            }
        }
        return std::pair<iterator, bool>(end(), false);
    }

    template<typename InputIterator>
    void insert(InputIterator begin, InputIterator end)
    {
        for (InputIterator it(begin); it != end; ++it) {
            insert(it->first, it->second);
        }
    }

    const_iterator find(const KEY &key) const
    {
        DIMN idx;
        if ((idx = key_to_index(key)) < dense_dimension()) {
            return const_iterator(*this, DENSE::find(idx));
        } else {
            return const_iterator(*this, SPARSE::find(key));
        }
    }

    iterator find(const KEY &key)
    {
        DIMN idx;
        if ((idx = key_to_index(key)) < dense_dimension()) {
            return iterator(*this, DENSE::find(idx));
        } else {
            return iterator(*this, SPARSE::find(key));
        }
    }

    void erase(const KEY &key)
    {
        DIMN idx;
        if ((idx = key_to_index(key)) < dense_dimension()) {
            DENSE::value(idx) = zero;
        } else {
            SPARSE::erase(key);
        }
    }

    void erase(iterator &it)
    {
        if (it->is_dense()) {
            it->value() = zero;
        } else {
            SPARSE::erase(it->key());
        }
    }


public:

    void clear()
    {
        DENSE::clear();
        SPARSE::clear();
    }

public:

    void swap(hybrid_vector &other)
    {
        DENSE::swap(other);
        SPARSE::swap(other);
    }

public:

    // Element access

    const SCALAR &operator[](const KEY &key) const
    {
        DIMN idx;
        if ((idx = key_to_index(key)) < DENSE::dimension()) {
            return DENSE::value(idx);
        } else {
            return SPARSE::operator[](key);
        }
    }

    SCALAR &operator[](const KEY &key)
    {
        DIMN idx;
        if ((idx = key_to_index(key)) < DENSE::dimension()) {
            return DENSE::value(idx);
        } else {
            return SPARSE::operator[](key);
        }
    }

    SCALAR &value(const DIMN idx)
    {
        if (idx < DENSE::dimension()) {
            return DENSE::value(idx);
        } else {
            return SPARSE::operator[](index_to_key(idx));
        }
    }

    const SCALAR &value(const DIMN idx) const
    {
        if (idx < DENSE::dimension()) {
            return DENSE::value(idx);
        } else {
            return SPARSE::operator[](index_to_key(idx));
        }
    }

    SCALAR &dense_value(const DIMN idx)
    {
        return DENSE::value(idx);
    }

    const SCALAR &dense_value(const DIMN idx) const
    {
        return DENSE::value(idx);
    }


public:

    // Comparison operators

    bool operator==(const hybrid_vector &rhs) const
    {
        if (empty() && rhs.empty()) {
            return true;
        }

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

        if (dense_part_eq.first == dense_dimension()
            && dense_part_eq.first == rhs.dense_dimension()) {
            return sparse_part() == rhs.sparse_part();
        } else if (dense_part_eq.first == dense_dimension()) {
            cend = SPARSE::end();

            for (DIMN i = dense_part_eq.first; i < rhs.dense_dimension(); ++i) {
                if (rhs.dense_part().value(i) == zero) {
                    continue;
                }

                cit = SPARSE::find(index_to_key(i));
                if (cit == cend || cit->value() != rhs.dense_part().value(i)) {
                    return false;
                }
            }

            for (oit = rhs.sparse_part().begin(); oit != rhs.sparse_part().end(); ++oit) {
                cit = SPARSE::find(oit->key());
                if (cit == cend || cit->value() != oit->value()) {
                    return false;
                }
            }

            return true;
        } else if (dense_part_eq.first == rhs.dense_dimension()) {

            cend = rhs.sparse_part().end();

            for (DIMN i = dense_part_eq.first; i < dense_dimension(); ++i) {
                if (DENSE::value(i) == zero) {
                    continue;
                }

                cit = rhs.sparse_part().find(index_to_key(i));
                if (cit == cend || cit->value() != DENSE::value(i)) {
                    return false;
                }
            }

            for (oit = sparse_part().begin(); oit != sparse_part().end(); ++oit) {
                cit = rhs.sparse_part().find(oit->key());
                if (cit == cend || cit->value() != oit->value()) {
                    return false;
                }
            }
            return true;

        }

        // Should be unreachable
        assert(false);
        return false;


    }

    bool operator!=(const hybrid_vector &other) const
    {
        return !operator==(other);
    }


public:

    // Arithmetic operators
    // The first few are operators that are simply applied to each
    // part individually.

    hybrid_vector operator-(void) const
    {
        return hybrid_vector(
                DENSE::operator-(),
                SPARSE::operator-()
        );
    }

    hybrid_vector &operator*=(const SCALAR s)
    {
        DENSE::operator*=(s);
        SPARSE::operator*=(s);
        return *this;
    }

    hybrid_vector &operator/=(const RATIONAL s)
    {
        DENSE::operator/=(s);
        SPARSE::operator/=(s);
        return *this;
    }


    // The next few operators require some matching

    hybrid_vector &operator+=(const hybrid_vector &other)
    {
        DIMN dim = std::max(
                dense_dimension(),
                other.dense_dimension()
        );
        resize_dense_to_dimension(dim);

        DENSE::operator+=(other);
        SPARSE::operator+=(other);
        incorporate_sparse();
        return *this;
    }

    hybrid_vector &operator-=(const hybrid_vector &other)
    {
        DIMN dim = std::max(
                dense_dimension(),
                other.dense_dimension()
        );
        std::cout << dim << '\n';
        resize_dense_to_dimension(dim);

        DENSE::operator-=(other);
        SPARSE::operator-=(other);
        incorporate_sparse();
        return *this;
    }

    hybrid_vector &operator&=(const hybrid_vector &rhs)
    {
        DIMN dim = std::max(
                dense_dimension(),
                rhs.dense_dimension()
        );
        resize_dense_to_dimension(dim);

        DENSE::operator&=(rhs);
        SPARSE::operator&=(rhs);
        incorporate_sparse();
        return *this;
    }

    hybrid_vector &operator|=(const hybrid_vector &rhs)
    {
        DIMN dim = std::max(
                dense_dimension(),
                rhs.dense_dimension()
        );
        resize_dense_to_dimension(dim);

        DENSE::operator|=(rhs);
        SPARSE::operator|=(rhs);
        incorporate_sparse();
        return *this;
    }


    __DECLARE_BINARY_OPERATOR(hybrid_vector, *, *=, SCALAR);

    __DECLARE_BINARY_OPERATOR(hybrid_vector, /, /=, RATIONAL);

    __DECLARE_BINARY_OPERATOR(hybrid_vector, +, +=, hybrid_vector);

    __DECLARE_BINARY_OPERATOR(hybrid_vector, -, -=, hybrid_vector);

    __DECLARE_BINARY_OPERATOR(hybrid_vector, &, &=, hybrid_vector);

    __DECLARE_BINARY_OPERATOR(hybrid_vector, |, |=, hybrid_vector);


public:

    // Fused operations

    DEFINE_FUSED_OP(add_scal_prod, SCALAR, +=, *);

    DEFINE_FUSED_OP(sub_scal_prod, SCALAR, -=, *);

    DEFINE_FUSED_OP(add_scal_div, RATIONAL, +=, /);

    DEFINE_FUSED_OP(sub_scal_div, RATIONAL, -=, /);


public:

    // Norms

    SCALAR NormL1() const
    {
        return DENSE::NormL1() + SPARSE::NormL1();
    }

    SCALAR NormL1(const DEG deg) const
    {
        return DENSE::NormL1(deg) + SPARSE::NormL1(deg);
    }

    SCALAR NormLInf() const
    {
        SCALAR dli = DENSE::NormLInf();
        SCALAR sli = SPARSE::NormLInf();
        return std::max(dli, sli);
    }

    SCALAR NormLInf(const DEG deg) const
    {
        return std::max(DENSE::NormLInf(deg), SPARSE::NormLInf(deg));
    }

public:

    // Display

    inline friend std::ostream &operator<<(std::ostream &os,
                                           const hybrid_vector &rhs)
    {
        std::pair<BASIS *, KEY> token;
        token.first = &basis;

        os << '{';

        for (const_iterator cit(rhs.begin()); cit != rhs.end(); ++cit) {
            if (zero != cit->value()) {
                token.second = cit->key();
                os << ' ' << cit->value() << '(' << token << ')';
            }
        }

        os << ' ' << '}';
        return os;
    }

public:

    static bool comp(std::pair<KEY, SCALAR> lhs,
                     std::pair<KEY, SCALAR> rhs)
    {
        key_ordering ord;
        return ord(lhs.first, rhs.first);
    }


private:

    template<typename Vector, typename KeyTransform>
    void triangular_binary_transform_mixed_cases(
            Vector &result,
            const hybrid_vector &rhs,
            KeyTransform key_transform,
            const DEG max_depth
    ) const
    {
        /*
         * At this stage there are three remaining cases:
         *   - dense * sparse;
         *   - sparse * dense;
         *   - sparse * sparse;
         * If both sparse parts are empty, there is nothing further to do,
         * so we can return. Otherwise, we should separate the rhs by degree
         * to avoid map lookups in a tight loop, and then proceed with the
         * multiplication as usual.
         */
        if (sparse_empty() && rhs.sparse_empty()) {
            return;
        }

        std::vector<std::pair<KEY, SCALAR> > buffer;
        std::vector<typename std::vector<std::pair<KEY, SCALAR> >::const_iterator>
                iterators;
        SPARSE::separate_by_degree(buffer, rhs, max_depth, iterators);

        // Do the dense * sparse first
        typename std::vector<std::pair<KEY, SCALAR> >::const_iterator cit,
                buf_begin(buffer.begin());

        DEG rh_max_deg;
        if (!dense_empty()) {
            for (DEG lhs_deg = 0; lhs_deg <= std::min(dense_degree(), max_depth); ++lhs_deg) {
                rh_max_deg = max_depth - lhs_deg;
                for (DIMN i = DENSE::start_of_degree(lhs_deg); i < DENSE::start_of_degree(lhs_deg + 1); ++i) {
                    for (cit = buf_begin; cit != iterators[rh_max_deg]; ++cit) {
                        key_transform(result, index_to_key(i), DENSE::value(i), cit->first, cit->second);
                    }
                }
            }
        }

        // Now do sparse * dense and sparse * sparse in one step.
        DEG lhs_deg;
        for (typename SPARSE::const_iterator it(SPARSE::begin()); it != SPARSE::end(); ++it) {
            lhs_deg = basis.degree(it->key());
            assert(lhs_deg <= max_depth);
            rh_max_deg = max_depth - lhs_deg;

            if (!rhs.dense_empty()) {
                for (DEG rhs_deg = 0; rhs_deg <= std::min(rh_max_deg, rhs.dense_degree()); ++rhs_deg) {
                    for (DIMN j = DENSE::start_of_degree(rhs_deg); j < DENSE::start_of_degree(rhs_deg + 1); ++j) {
                        key_transform(result, it->key(), it->value(), index_to_key(j), rhs.dense_value(j));
                    }
                }
            }

            for (cit = buf_begin; cit != iterators[rh_max_deg]; ++cit) {
                key_transform(result, it->key(), it->value(), cit->first, cit->second);
            }

        }

        dtl::vector_base_access::convert(result).incorporate_sparse();
        //dtl::vector_base_access::convert(result).maybe_resize();
    }

    template<typename Vector, typename KeyTransform>
    void square_binary_transform_mixed_case(
            Vector &result,
            const hybrid_vector &rhs,
            KeyTransform key_transform
    ) const
    {
        /*
         * At this stage there are three remaining cases:
         *   - dense * sparse;
         *   - sparse * dense;
         *   - sparse * sparse;
         * If both sparse parts are empty, there is nothing further to do,
         * so we can return. Otherwise, we should separate the rhs by degree
         * to avoid map lookups in a tight loop, and then proceed with the
         * multiplication as usual.
         */
        if (sparse_empty() && rhs.sparse_empty()) {
            return;
        }
        std::vector<std::pair<KEY, SCALAR> > buffer;
        rhs.fill_buffer(buffer);

        typename std::vector<std::pair<KEY, SCALAR> >::const_iterator
                cit, buf_begin(buffer.begin()), buf_end(buffer.end());

        // First deal with the dense * sparse case
        for (DIMN i = 0; i < dense_dimension(); ++i) {
            for (cit = buf_begin; cit != buf_end; ++cit) {
                key_transform(result, index_to_key(i), dense_value(i), cit->first, cit->second);
            }
        }

        // Now deal with sparse * dense and sparse * sparse together
        for (typename SPARSE::const_iterator it(SPARSE::begin()); it != SPARSE::end(); ++it) {
            for (DIMN j = 0; j < rhs.dense_dimension(); ++j) {
                key_transform(result, it->key(), it->value(), index_to_key(j), rhs.dense_value(j));
            }

            for (cit = buf_begin; cit != buf_end; ++cit) {
                key_transform(result, it->key(), it->value(), cit->first, cit->second);
            }
        }

        dtl::vector_base_access::convert(result).maybe_resize();
    }

public:

    // Transform methods

    template<typename Vector, typename KeyTransform>
    void triangular_buffered_apply_binary_transform(
            Vector &result,
            const hybrid_vector &rhs,
            KeyTransform key_transform,
            const DEG max_depth
    ) const
    {

        if (!dense_empty() && !rhs.dense_empty()) {
            DEG max_degree = std::min(max_depth, dense_degree() + rhs.dense_degree());
            DENSE::triangular_buffered_apply_binary_transform(
                    result, rhs, key_transform, max_degree
            );
        }

        triangular_binary_transform_mixed_cases(result, rhs, key_transform, max_depth);

    }

    template<typename Vector, typename KeyTransform, typename IndexTransform>
    void triangular_buffered_apply_binary_transform(
            Vector &result,
            const hybrid_vector &rhs,
            KeyTransform key_transform,
            IndexTransform index_transform,
            const DEG max_depth
    ) const
    {
        // The dense part will be resized to accommodate the max_dense_degree
        DEG max_dense_degree = std::min(max_depth, dense_degree() + rhs.dense_degree());

        if (!dense_empty() && !rhs.dense_empty()) {
            DENSE::triangular_buffered_apply_binary_transform(
                    result, rhs, key_transform, index_transform, max_dense_degree
            );
        }

        triangular_binary_transform_mixed_cases(result, rhs, key_transform, max_depth);

    }

    template<typename Vector, typename KeyTransform>
    void square_buffered_apply_binary_transform(
            Vector &result,
            const hybrid_vector &rhs,
            KeyTransform key_transform
    ) const
    {
        if (!dense_empty() && !rhs.dense_empty()) {
            DENSE::square_buffered_apply_binary_transform(result, rhs, key_transform);
        }

        square_binary_transform_mixed_case(result, rhs, key_transform);
    }

    template<typename Vector, typename KeyTransform, typename IndexTransform>
    void square_buffered_apply_binary_transform(
            Vector &result,
            const hybrid_vector &rhs,
            KeyTransform key_transform,
            IndexTransform index_transform
    ) const
    {
        if (!dense_empty() && !rhs.dense_empty()) {
            DENSE::square_buffered_apply_binary_transform(result, rhs, key_transform, index_transform);
        }

        square_binary_transform_mixed_case(result, rhs, key_transform);
    }



public:

    template <typename Transform>
    void buffered_apply_unary_transform(
            hybrid_vector& result,
            Transform transform,
            const DEG max_deg
    ) const
    {
        if (empty()) {
            return;
        }

        if (!dense_empty()) {
            DENSE::buffered_apply_unary_transform(result, transform, max_deg);
        }

        SPARSE::buffered_apply_unary_transform(result, transform, max_deg);
        result.maybe_resize();
    }


    template <typename Transform>
    void buffered_apply_unary_transform(
            hybrid_vector& result,
            Transform transform
    ) const
    {
        if (empty()) {
            return;
        }

        if (!dense_empty()) {
            DENSE::buffered_apply_unary_transform(result, transform);
        }
        SPARSE::buffered_apply_unary_transform(result, transform);
        result.maybe_resize();
    }





};

#undef DEFINE_FUSED_OP


} // namespace vectors
} // namespace alg





#endif //LIBALGEBRA_HYBRID_VECTOR_H
