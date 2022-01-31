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

    using pair_type = std::pair<key_type, scalar_type>;
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

private:
    explicit flat_sparse_vector(storage_type&& storage) : m_storage(std::move(storage))
    {}

    pair_type* find_key(const key_type& key)
    {
        return std::binary_search(m_storage.begin(), m_storage.end(),
                                  pair_type(key, zero),
                                  typename basis_type::ordering_tag::pair_order());
    }

    const pair_type* find_key(const key_type& key) const
    {
        return std::binary_search(m_storage.begin(), m_storage.end(),
                                  pair_type(key, zero),
                                  typename basis_type::ordering_tag::pair_order());
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
    // Swap operation

    void swap(flat_sparse_vector& rhs)
    {
        std::swap(m_storage, rhs.m_storage);
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
    storage_type apply_operation(const flat_sparse_vector& rhs, Fn&& function)
    {
        const pair_type *lhs_ptr = m_storage.begin(), lhs_end = m_storage.end();
        const pair_type *rhs_ptr = rhs.m_storage.begin(), rhs_end = rhs.m_storage.end();

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
            for (auto item : m_storage) {
                if (item->first != key) {
                    new_storage.emplace(idx++, item);
                }
            }
            m_storage = std::move(new_storage);
        }
        else {
            storage_type new_storage;
            new_storage.reserve(m_storage.size() + 1);
            typename basis_type::ordering_tag::order order;

            size_type idx = 0;
            auto it = m_storage.begin();
            for (; it != m_storage.end() && order(it->first, key); ++it) {
                new_storage.emplace(idx++, *it);
            }

            new_storage.emplace(idx++, {key, function(zero)});

            for (; it != m_storage.end(); ++it) {
                new_storage.emplace(idx++, *it);
            }
            m_storage = std::move(new_storage);
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
        for (auto item : rhs.m_storage) {
            token.second = item->first;
            os << ' ' << item->second << '(' << token << ')';
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

    DEG degree() const
    {
        DEG ans(0);
        for (auto item : m_storage) {
            ans = std::max(ans, base_vector_type::basis.degree(item->first));
        }
        return ans;
    }

    bool degree_equals(const DEG degree) const
    {
        bool result(false);
        for (auto item : m_storage) {
            auto deg = base_vector_type::basis.degree(item->first);
            if (deg > degree) {
                return false;
            } else if (!result && deg == degree) {
                result = true;
            }
        }
        return result;
    }


};

}// namespace vectors
}// namespace alg

#endif//LIBALGEBRA_FLAT_SPARSE_VECTOR_H
