//
// Created by sam on 23/09/2021.
//

#ifndef LIBALGEBRA_HALL_SET_H
#define LIBALGEBRA_HALL_SET_H

#include <cassert>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

#include "implementation_types.h"

namespace alg {

template <DEG Width>
class hall_set
{

    using degree_type               = DEG;
    using letter_type               = LET;
    using key_type                  = letter_type;
    using size_type                 = DIMN;

    using parent_type               = std::pair<key_type, key_type>;
    using reverse_map_type          = std::map<parent_type, key_type>;
    using data_type                 = std::vector<parent_type>;
    using letter_vec_type           = std::vector<letter_type>;
    using l2k_map_type              = std::vector<key_type>;

    letter_vec_type             m_letters;
    data_type                   m_data;
    reverse_map_type            m_reverse_map;
    degree_type                 m_current_degree;
    l2k_map_type                m_l2k;





private:

    hall_set()
    {

    }

public:

    static hall_set& get_instance(degree_type degree)
    {
        static hall_set data;
        static std::mutex lock;
        std::lock_guard<std::mutex> access(lock);

        if (data.m_current_degree < degree) {
            data.grow_up(degree);
        }

        return data;
    }


private:

    void grow_up(degree_type degree)
    {
        for (degree_type d = m_current_degree; d <= degree; ++d) {
            for (degree_type e=1; 2*e <= d; ++e) {

            }

            ++m_current_degree;
        }
    }


public:


    key_type begin() const noexcept
    {
        return key_type(1);
    }

    key_type end(degree_type width) const noexcept
    {
        return key_type(0);
    }

    key_type next_key(const key_type& key) const noexcept
    {
        if (key < size()) {
            return (key+1);
        } else {
            return key_type(0);
        }
    }


private:

    static std::string letter_to_string(letter_type letter) noexcept
    {
        return std::to_string(letter);
    }

    static std::string key2string_binop(const std::string& a, const std::string& b)
    {
        return "[" + a + "," + b + "]";
    }


public:


    key_type key_of_letter(letter_type letter) const
    {
        assert(letter(key));
        return m_l2k[letter-1];
    }

    key_type lparent(const key_type& key) const noexcept
    {
        return m_data[key].first;
    }

    key_type rparent(const key_type& key) const noexcept
    {
        return m_data[key].second;
    }

    bool letter(const key_type& key) const noexcept
    {
        return ((key > 0) && (key <= Width));
    }

    size_type size() const noexcept
    {
        return m_data.size() - 1;
    }

public:

    struct no_caching_tag { };

    template<DEG CacheDepth>
    struct lazy_cache_tag { };

    template<DEG CacheDepth>
    struct lookup_table_tag { };

    template<typename Function, typename BinOp, typename Tag>
    class extended_function {
    public:
        typedef Function function_type;
        typedef BinOp binary_operation_type;
        typedef Tag tag_type;
        typedef decltype(std::declval<function_type>()(std::declval<key_type>())) output_type;

        using table_t = std::unordered_map<key_type, output_type>;

        extended_function()
                : m_hall_set(hall_set::get_instance(0)), m_tag(), m_fn(), m_op() { }

        extended_function(Function fn, BinOp op)
                : m_hall_set(hall_set::get_instance(0)), m_tag(), m_fn(fn), m_op(op) { }

    private:

        output_type eval_impl(const key_type& k) const
        {
            if (m_hall_set.letter(k)) {
                return m_fn(m_hall_set.getletter(k));
            }
            else {
                return m_op(
                        operator()(m_hall_set.lparent(k)),
                        operator()(m_hall_set.rparent(k))
                );
            }
        }

        output_type eval(const key_type& k, no_caching_tag) const
        {
            return eval_impl(k);
        }

        template<DEG CacheDepth>
        output_type eval(const key_type& k, lazy_cache_tag<CacheDepth>) const
        {
            static boost::recursive_mutex table_lock;
            static table_t table;

            if (m_hall_set.degree(k)<=CacheDepth) {
                boost::lock_guard<boost::recursive_mutex> access(table_lock);

                typename table_t::iterator it = table.find(k);
                if (it!=table.end()) {
                    return it->second;
                }

                return table[k] = eval_impl(k);
            }
            else {
                return eval_impl(k);
            }
        }

        output_type eval(const key_type& k, lazy_cache_tag<0>) const
        {
            static boost::recursive_mutex table_lock;
            static table_t table;

            boost::lock_guard<boost::recursive_mutex> access(table_lock);

            typename table_t::iterator it = table.find(k);
            if (it!=table.end()) {
                return it->second;
            }

            return table[k] = eval_impl(k);
        }

        template<DEG Depth>
        table_t fill_table() const
        {
            table_t result;

            key_type k = 1;
            for (; m_hall_set.degree(k)==1; ++k) {
                result[k] = m_fn(m_hall_set.getletter(k));
            }

            for (; m_hall_set.degree(k)<=Depth; ++k) {
                result[k] = m_op(result[m_hall_set.lparent(k)], result[m_hall_set.rparent(k)]);
            }

            return result;
        }

        template<DEG CacheDepth>
        output_type eval(const key_type& k, lookup_table_tag<CacheDepth>) const
        {
            static table_t table = fill_table<CacheDepth>();

            if (m_hall_set.degree(k)<=CacheDepth) {
                return table[k];
            }
            else {
                return eval_impl(k);
            }
        }

    public:

        output_type operator()(const key_type& k) const
        {
            return eval(k, m_tag);
        }

    private:
        hall_set& m_hall_set;
        function_type m_fn;
        binary_operation_type m_op;
        tag_type m_tag;
    };


    using key2string_type = extended_function<
            decltype(&hall_set::letter_to_string),
            decltype(&hall_set::key2string_binop),
            lookup_table_tag<0>
    >;

    static const key2string_type key2string;

};

template <DEG Width>
const typename hall_set<Width>::key2string_type hall_set<Width>::key2string(
        hall_set<Width>::letter_to_string, hall_set<Width>::key2string_binop
        );






}




#endif //LIBALGEBRA_HALL_SET_H
