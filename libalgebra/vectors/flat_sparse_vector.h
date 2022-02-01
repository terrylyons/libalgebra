//
// Created by sam on 31/01/2022.
//

#ifndef LIBALGEBRA_FLAT_SPARSE_VECTOR_H
#define LIBALGEBRA_FLAT_SPARSE_VECTOR_H

#include <libalgebra/implementation_types.h>
#include <libalgebra/utils/order_trait.h>
#include <libalgebra/vectors/base_vector.h>
#include <libalgebra/vectors/dense_storage.h>
#include <libalgebra/vectors/iterators.h>

#include <algorithm>
#include <iosfwd>
#include <map>
#include <utility>
#include <vector>

namespace alg {
namespace vectors {

template<typename Basis, typename Coeffs>
class flat_sparse_vector : public base_vector<Basis, Coeffs>
{
    using basis_type = Basis;
    using key_type = typename Basis::KEY;
    using scalar_type = typename Coeffs::SCA;

    using pair_type = std::pair<const key_type, scalar_type>;
    using storage_type = dense_storage<pair_type>;

    using base_vector_type = base_vector<Basis, Coeffs>;
    using coeff_type = Coeffs;

    class iterator_item;
    class const_iterator_item;

    using size_type = typename storage_type::size_type;

public:
    // Typedefs for compatibility with the library
    using BASIS = Basis;
    using KEY = key_type;
    using SCALAR = scalar_type;
    using RATIONAL = typename Coeffs::Q;

    using base_vector_type::zero;

    using iterator = iterators::vector_iterator<iterator_item>;
    using const_iterator = iterators::vector_iterator<const_iterator_item>;

    using key_ordering = typename dtl::requires_order<Basis>::key_ordering;

private:
    storage_type m_storage;

public:
    // Constructors - We're only defining non-trivial constructors that won't
    // be filled in by the compiler

    flat_sparse_vector() : m_storage()
    {}

    explicit flat_sparse_vector(const key_type& key, const scalar_type& scalar = base_vector_type::one)
        : m_storage{{key, scalar}}
    {}

    explicit flat_sparse_vector(const scalar_type* it, const scalar_type* end)
        : m_storage(static_cast<size_type>(end - it))
    {
        auto key = base_vector_type::basis.begin();
        auto key_end = base_vector_type::basis.end();
        size_type idx = 0;
        while (it != end && key != key_end) {
            m_storage.emplace(idx++, pair_type(key, *(it++)));
            key = base_vector_type::basis.nextkey(key);
        }
    }

    explicit flat_sparse_vector(size_type offset, const scalar_type* it, const scalar_type* end)
        : m_storage(static_cast<size_type>(end - it))
    {
        auto key = base_vector_type::basis.begin();
        auto key_end = base_vector_type::basis.end();

        for (size_type i = 0; i < offset; ++i) {
            key = base_vector_type::basis.nextkey(key);
        }

        size_type idx = 0;
        while (it != end && key != key_end) {
            m_storage.emplace(idx++, pair_type(key, *(it++)));
            key = base_vector_type::basis.nextkey(key);
        }
    }

public:
    // Vector information methods

    size_type size() const
    {
        return m_storage.size();
    }

    DEG degree() const
    {
        if (m_storage.empty()) {
            return 0;
        }
        // Short cut assuming that the ordering respects degree. Will have to be changed
        return base_vector_type::basis.degree(m_storage[m_storage.size() - 1].first);
    }

    /*
    DEG degree() const
    {
        DEG ans(0);
        for (auto item : m_storage) {
            ans = std::max(ans, base_vector_type::basis.degree(item->first));
        }
        return ans;
    }
     */

    bool empty() const
    {
        return m_storage.empty();
    }

private:
    explicit flat_sparse_vector(storage_type&& storage) : m_storage(std::move(storage))
    {}

    pair_type* find_key(const key_type& key)
    {
        return std::find_if(m_storage.begin(), m_storage.end(),
                            [&key](pair_type& arg) {
                                return arg.first == key;
                            });
    }

    const pair_type* find_key(const key_type& key) const
    {
        return std::find_if(m_storage.begin(), m_storage.end(),
                            [&key](const pair_type& arg) {
                                return arg.first == key;
                            });
    }

public:
    // Element access

    scalar_type& operator[](const key_type& key)
    {
        auto* found = find_key(key);
        if (found != m_storage.end()) {
            return found->second;
        }
        else {
            throw std::invalid_argument("no value corresponding to key");
        }
    }

    const scalar_type& operator[](const key_type& key) const
    {
        auto* found = find_key(key);
        if (found != m_storage.end()) {
            return found->second;
        }
        else {
            return zero;
        }
    }

public:
    // Iterator interaction methods

    iterator begin()
    {
        return iterator(*this, m_storage.begin());
    }

    iterator end()
    {
        return iterator(*this, m_storage.end());
    }

    const_iterator begin() const
    {
        return const_iterator(*this, m_storage.begin());
    }

    const_iterator end() const
    {
        return const_iterator(*this, m_storage.end());
    }

    const_iterator cbegin() const
    {
        return begin();
    }

    const_iterator cend() const
    {
        return end();
    }

    iterator find(const key_type& key)
    {
        return iterator(*this, find_key(key));
    }

    const_iterator find(const key_type& key) const
    {
        return const_iterator(*this, find_key(key));
    }

private:
    template<typename SortedInputIt>
    storage_type copy_insert(size_type count, SortedInputIt cit, SortedInputIt cend) const
    {
        storage_type new_storage;
        new_storage.reserve(m_storage.size() + count);
        typename basis_type::ordering_tag::order order;

        size_type idx = 0;
        auto it = m_storage.begin();
        for (; cit != cend; ++cit) {
            for (; it != m_storage.end() && order(it->first, cit->first); ++it) {
                new_storage.emplace(idx++, *it);
            }

            new_storage.emplace(idx++, {cit->first, cit->second});
        }

        for (; it != m_storage.end(); ++it) {
            new_storage.emplace(idx++, *it);
        }
        return new_storage;
    }

    storage_type copy_insert(std::initializer_list<pair_type> args) const
    {
        return copy_insert(args.size(), args.begin(), args.end());
    }

public:
    void insert(const key_type& key, scalar_type val)
    {
        auto* found = find_key(key);
        if (found != m_storage.end()) {
            found->second = std::move(val);
        }
        else {
            m_storage = copy_insert({{key, val}});
        }
    }

    std::pair<iterator, bool> insert(pair_type& arg)
    {
        auto* found = find_key(arg.first);
        if (found != m_storage.end()) {
            return {iterator(*this, found), false};
        }
        else {
            m_storage = copy_insert({arg});
        }
    }

    template<typename InputIt>
    void insert(InputIt begin, InputIt end)
    {
        std::vector<pair_type> buffer(begin, end);
        std::sort(buffer.begin(), buffer.end(), typename basis_type::ordering_type::pair_order());

        m_storage = copy_insert(buffer.size(), buffer.begin(), buffer.end());
    }

private:
    template<typename SortedInputIt>
    storage_type copy_erase(size_type count, SortedInputIt cit, SortedInputIt cend) const
    {
        storage_type new_storage;
        new_storage.reserve(m_storage.size() + count);
        typename basis_type::ordering_tag::order order;

        size_type idx = 0;
        auto it = m_storage.begin();
        for (; cit != cend; ++cit) {
            for (; it != m_storage.end() && order(it->first, *cit); ++it) {
                new_storage.emplace(idx++, *it);
            }
            ++it;
        }

        return new_storage;
    }

    storage_type copy_erase(std::initializer_list<key_type> args) const
    {
        std::sort(args.begin(), args.end(), typename basis_type::ordering_tag::order());
        copy_erase(args.size(), args.begin(), args.end());
    }

public:
    void erase(iterator& it)
    {
        m_storage = copy_erase({it->first});
    }

    void erase(const key_type& key)
    {
        m_storage = copy_erase({key});
    }

    void clear()
    {
        m_storage.clear();
    }

public:
    // Swap operation

    void swap(flat_sparse_vector& rhs)
    {
        std::swap(m_storage, rhs.m_storage);
    }

public:
    // Unary minus

    flat_sparse_vector operator-() const
    {
        storage_type new_storage;
        new_storage.reserve(m_storage.size());

        size_type idx = 0;
        for (auto& item : m_storage) {
            new_storage.emplace(idx++, {item.first, coeff_type::uminus(item.second)});
        }

        return flat_sparse_vector(std::move(new_storage));
    }

public:
    // Inplace arithmetic

    flat_sparse_vector& operator*=(const scalar_type& scalar)
    {
        if (scalar == zero) {
            m_storage.resize(0);
        }
        else {
            for (auto& pair : m_storage) {
                coeff_type::mul_inplace(pair.second, scalar);
            }
        }
        return *this;
    }

    flat_sparse_vector& operator/=(const RATIONAL& rational)
    {
        for (auto& pair : m_storage) {
            coeff_type::div_inplace(pair.second, rational);
        }
        return *this;
    }

private:
    template<typename Fn>
    storage_type apply_operation(const flat_sparse_vector& rhs, Fn&& function) const
    {
        const pair_type *lhs_ptr = m_storage.begin(), *lhs_end = m_storage.end();
        const pair_type *rhs_ptr = rhs.m_storage.begin(), *rhs_end = rhs.m_storage.end();

        std::vector<pair_type> buffer;
        buffer.reserve(m_storage.size() + rhs.m_storage.size());
        while (lhs_ptr != lhs_end && rhs_ptr != rhs_end) {
            if (lhs_ptr->first == rhs_ptr->first) {
                auto result = function(lhs_ptr->second, rhs_ptr->second);
                if (result != zero) {
                    buffer.emplace_back(lhs_ptr->first, result);
                }
                ++lhs_ptr;
                ++rhs_ptr;
            }
            else if (lhs_ptr->first < rhs_ptr->first) {
                buffer.emplace_back(lhs_ptr->first, function(lhs_ptr->second, zero));
                ++lhs_ptr;
            }
            else {
                buffer.emplace_back(rhs_ptr->first, function(zero, rhs_ptr->second));
                ++rhs_ptr;
            }
        }

        for (; lhs_ptr != lhs_end; ++lhs_ptr) {
            buffer.emplace_back(lhs_ptr->first, function(lhs_ptr->second, zero));
        }

        for (; rhs_ptr != rhs_end; ++rhs_ptr) {
            buffer.emplace_back(rhs_ptr->first, function(zero, rhs_ptr->second));
        }

        return storage_type(buffer.data(), buffer.data() + buffer.size());
    }

public:
    flat_sparse_vector& operator+=(const flat_sparse_vector& rhs)
    {
        m_storage = apply_operation(rhs, coeff_type::add);
        return *this;
    }

    flat_sparse_vector& operator-=(const flat_sparse_vector& rhs)
    {
        m_storage = apply_operation(rhs, coeff_type::sub);
        return *this;
    }

    // External binary operations

    flat_sparse_vector operator*(const scalar_type& scalar) const
    {
        if (scalar == zero) {
            return flat_sparse_vector();
        }

        storage_type result_storage;
        result_storage.reserve(m_storage.size());

        size_type idx = 0;
        for (auto item : m_storage) {
            result_storage.emplace(idx++, {item.first, coeff_type::mul(item.second, scalar)});
        }

        return flat_sparse_vector(result_storage);
    }

    flat_sparse_vector operator/(const RATIONAL& rational) const
    {
        storage_type result_storage;
        result_storage.reserve(m_storage.size());

        size_type idx = 0;
        for (auto item : m_storage) {
            result_storage.emplace(idx++, {item.first, coeff_type::div(item.second, rational)});
        }

        return flat_sparse_vector(result_storage);
    }

public:
    flat_sparse_vector operator+(const flat_sparse_vector& rhs) const
    {
        return flat_sparse_vector(apply_operation(rhs, coeff_type::add));
    }

    flat_sparse_vector operator-(const flat_sparse_vector& rhs) const
    {
        return flat_sparse_Vector(apply_operation(rhs, coeff_type::sub));
    }

    // Fused operations

    flat_sparse_vector& add_scal_prod(const flat_sparse_vector& rhs, const scalar_type& scalar)
    {
        if (scalar == zero || rhs.m_storage.empty()) {
            return *this;
        }

        m_storage = apply_operation(rhs, [&scalar](const scalar_type& left, const scalar_type& right) {
            return coeff_type::add(left, coeff_type::mul(scalar, right));
        });

        return *this;
    }

    flat_sparse_vector& sub_scal_prod(const flat_sparse_vector& rhs, const scalar_type& scalar)
    {
        if (scalar == zero || rhs.m_storage.empty()) {
            return *this;
        }

        m_storage = apply_operation(rhs, [&scalar](const scalar_type& left, const scalar_type& right) {
            return coeff_type::sub(left, coeff_type::mul(scalar, right));
        });

        return *this;
    }

    flat_sparse_vector& add_scal_div(const flat_sparse_vector& rhs, const RATIONAL& rational)
    {
        if (rhs.m_storage.empty()) {
            return *this;
        }

        m_storage = apply_operation(rhs, [&rational](const scalar_type& left, const scalar_type& right) {
            return coeff_type::add(left, coeff_type::div(right, rational));
        });

        return *this;
    }

    flat_sparse_vector& sub_scal_div(const flat_sparse_vector& rhs, const RATIONAL& rational)
    {
        if (rhs.m_storage.empty()) {
            return *this;
        }

        m_storage = apply_operation(rhs, [&rational](const scalar_type& left, const scalar_type& right) {
            return coeff_type::sub(left, coeff_type::div(right, rational));
        });

        return *this;
    }

private:
    // fused operations using keys implementation

    template<typename Fn>
    void apply_fused_key_op(const key_type& key, Fn function)
    {
        auto* found = find_key(key);
        if (found != m_storage.end() && (found->second = function(found->second)) != zero) {
            return;
        }
        else if (found != m_storage.end()) {
            storage_type new_storage;
            new_storage.reserve(m_storage.size());
            size_type idx = 0;
            for (auto& item : m_storage) {
                if (item.first != key) {
                    new_storage.emplace(idx++, item);
                }
            }
            m_storage = std::move(new_storage);
        }
        else {
            m_storage = copy_insert({{key, function(zero)}});
        }
    }

public:
    flat_sparse_vector& add_scal_prod(const key_type& rhs, const scalar_type& scalar)
    {
        apply_fused_key_op(rhs, [&scalar](const scalar_type& current) {
            return coeff_type::add(current, scalar);
        });
        return *this;
    }

    flat_sparse_vector& add_scal_div(const key_type& rhs, const RATIONAL& rational)
    {
        apply_fused_key_op(rhs, [&rational](const scalar_type& current) {
            return coeff_type::add(current, coeff_type::div(coeff_type::one, rational));
        });
        return *this;
    }

    flat_sparse_vector& sub_scal_prod(const key_type& rhs, const scalar_type& scalar)
    {
        apply_fused_key_op(rhs, [&scalar](const scalar_type& current) {
            return coeff_type::sub(current, scalar);
        });
        return *this;
    }

    flat_sparse_vector& sub_scal_div(const key_type& rhs, const RATIONAL& rational)
    {
        apply_fused_key_op(rhs, [&rational](const scalar_type& current) {
            return coeff_type::sub(current, coeff_type::div(coeff_type::one, rational));
        });
        return *this;
    }

public:
    // comparison operators

    bool operator==(const flat_sparse_vector& rhs) const
    {
        auto *lhs_ptr = m_storage.begin(), lhs_end = m_storage.end();
        auto *rhs_ptr = rhs.m_storage.begin(), rhs_end = rhs.m_storage.end();

        while (lhs_ptr != lhs_end && rhs_ptr != rhs_end) {
            if (lhs_ptr->first == rhs_ptr->first) {
                if (lhs_ptr->second != rhs_ptr->second) {
                    return false;
                }
            }
            else {
                return false;
            }
            lhs_ptr++;
            rhs_ptr++;
        }

        return lhs_ptr == lhs_end && rhs_ptr == rhs_end;
    }

    bool operator!=(const flat_sparse_vector& rhs) const
    {
        return !operator==(rhs);
    }

public:
    /// Stream out operator
    friend std::ostream& operator<<(std::ostream& os, const flat_sparse_vector& rhs)
    {
        std::pair<basis_type*, key_type> token;
        token.first = &base_vector_type::basis;

        os << '{';
        for (auto& item : rhs.m_storage) {
            token.second = item.first;
            os << ' ' << item.second << '(' << token << ')';
        }

        return os << " }";
    }

public:
    // Norms

    scalar_type NormL1() const
    {
        scalar_type ans(zero);
        for (auto item : m_storage) {
            coeff_type::add_inplace(ans, abs(item->second));
        }
        return abs;
    }

    scalar_type NormL1(DEG target_deg) const
    {
        scalar_type ans(zero);
        for (auto item : m_storage) {
            if (base_vector_type::basis.degree(item->first) == target_deg) {
                coeff_type::add_inplace(ans, abs(item->second));
            }
        }
        return ans;
    }

    scalar_type NormLInf() const
    {
        scalar_type ans(zero);
        for (auto item : m_storage) {
            ans = std::max(ans, abs(item->second));
        }
        return abs;
    }

    scalar_type NormLInf(DEG target_deg) const
    {
        scalar_type ans(zero);
        for (auto item : m_storage) {
            if (base_vector_type::basis.degree(item->first) == target_deg) {
                ans = std::max(ans, abs(item->second));
            }
        }
        return ans;
    }

public:
    bool degree_equals(const DEG degree) const
    {
        bool result(false);
        for (auto item : m_storage) {
            auto deg = base_vector_type::basis.degree(item->first);
            if (deg > degree) {
                return false;
            }
            else if (!result && deg == degree) {
                result = true;
            }
        }
        return result;
    }
};

template<typename Basis, typename Coeffs>
class flat_sparse_vector<Basis, Coeffs>::iterator_item
{
    friend class iterators::vector_iterator<iterator_item>;
    friend class flat_sparse_vector;

private:
    pair_type* m_ptr;

public:
    iterator_item(flat_sparse_vector& vect, pair_type* ptr) : m_ptr(ptr)
    {}

    key_type key()
    {
        return m_ptr->first;
    }

    scalar_type& value()
    {
        return m_ptr->second;
    }

private:
    bool compare_iterators(const iterator_item& other) const
    {
        return (m_ptr == other.m_ptr);
    }

    void advance()
    {
        ++m_ptr;
    }

public:
    DIMN index() const
    {
        return base_vector_type::basis.key_to_index(m_ptr->first);
    }
};

template<typename Basis, typename Coeffs>
class flat_sparse_vector<Basis, Coeffs>::const_iterator_item
{
    friend class iterators::vector_iterator<const_iterator_item>;
    friend class flat_sparse_vector;

private:
    const pair_type* m_ptr;

public:
    const_iterator_item(const flat_sparse_vector& vect, const pair_type* ptr) : m_ptr(ptr)
    {}

    key_type key()
    {
        return m_ptr->first;
    }

    const scalar_type& value()
    {
        return m_ptr->second;
    }

private:
    bool compare_iterators(const const_iterator_item& other) const
    {
        return (m_ptr == other.m_ptr);
    }

    void advance()
    {
        ++m_ptr;
    }

public:
    DIMN index() const
    {
        return base_vector_type::basis.key_to_index(m_ptr->first);
    }
};

}// namespace vectors
}// namespace alg

#endif//LIBALGEBRA_FLAT_SPARSE_VECTOR_H
