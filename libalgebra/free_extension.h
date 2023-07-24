//
// Created by sam on 22/11/2021.
//

#ifndef LIBALGEBRA_FREE_EXTENSION_H
#define LIBALGEBRA_FREE_EXTENSION_H

#include <map>
#include <mutex>
#include <type_traits>
#include <utility>

#include "detail/caching_tags.h"

namespace alg {
namespace operators {

template<typename Tag, typename Basis, typename Function>
class free_extension_operator_impl
{
public:
    using basis_type = Basis;
    using tag_type = Tag;
    using function_type = Function;
    using result_type = decltype(std::declval<Function>()(std::declval<typename Basis::KEY>()));
    using key_type = typename basis_type::KEY;

    using table_t = std::map<key_type, result_type>;

private:
    template<typename Predicate>
    result_type eval(const key_type& k, lazy_cache_tag<Predicate> tag) const
    {
        static std::recursive_mutex table_lock;
        static table_t table;

        if (tag.predicate(k)) {
            std::lock_guard<std::recursive_mutex> access(table_lock);

            typename table_t::iterator it = table.find(k);
            if (it != table.end()) {
                return it->second;
            }

            return table[k] = m_fn(k);
        }
        else {
            return m_fn(k);
        }
    }

    result_type eval(const key_type& k, lazy_cache_tag<void> tag) const
    {
        static std::recursive_mutex table_lock;
        static table_t table;
        std::lock_guard<std::recursive_mutex> access(table_lock);

        typename table_t::iterator it = table.find(k);
        if (it != table.end()) {
            return it->second;
        }

        return table[k] = m_fn(k);
    }

    template<typename Predicate>
    table_t fill_table(Predicate predicate) const
    {
        table_t result;
        basis_type basis;

        for (key_type k = basis.begin(); k != basis.end() && predicate(k); k = basis.nextkey(k)) {
            result[k] = m_fn(k);
        }

        return result;
    }

    template<typename Predicate>
    result_type eval(const key_type& k, lookup_table_tag<Predicate>) const
    {
        static const table_t table = fill_table(m_tag.predicate);

        if (m_tag.predicate(k)) {
            return table[k];
        }
        else {
            return m_fn(k);
        }
    }

    result_type eval(const key_type& k, no_caching_tag) const
    {
        return m_fn(k);
    }

public:
    template<typename Coeff, template<typename, typename, typename...> class VT>
    result_type operator()(const vectors::vector<Basis, Coeff, VT>& arg) const
    {
        result_type result;
        for (auto item : arg) {
            result.add_scal_prod(eval(item.key(), m_tag), item.value());
        }
        return result;
    }

private:
    tag_type m_tag;
    function_type m_fn;
};

}// namespace operators
}// namespace alg

#endif//LIBALGEBRA_FREE_EXTENSION_H
