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
#include <map>

#include "implementation_types.h"
#include "utils/integer_maths.h"

namespace alg {

// Tag types for extension functions caching

/**
 * @brief No caching tag
 *
 * Using this tag indicates that no caching should be used
 * and all values should be computed recursively to leaves
 * in all computations. This is most suitable for functions
 * that are computed rarely and are not expensive.
 */
struct no_caching_tag { };


/**
 * @brief Lazy caching tag
 *
 * Using this tag indicates that values should be cached on
 * first computation conditionally based on a predicate.
 * The predicate should be a function object that evaluates
 * true on keys that should be added to the cache and false
 * if the value should be computed dynamically every time.
 *
 * Passing a void predicate will enable lazy caching for
 * all values.
 *
 * @tparam Predicate Predicate function object
 */
template<typename Predicate>
struct lazy_cache_tag {
    Predicate predicate;
};

template<>
struct lazy_cache_tag<void> { };


/**
 * @brief Lookup table caching tag
 *
 * Using this tag indicates that values should be computed
 * and stored in a lookup table on first call. The tag takes a
 * predicate function object, like the lazy caching tag, which
 * determines the keys that should have their values written in
 * the table. Note that the table construction will terminate at
 * the first key which fails the predicate.
 *
 * @tparam Predicate Predicate function object
 */
template<typename Predicate>
struct lookup_table_tag {
    Predicate predicate;
};




/**
 * @brief Data content of a Hall set
 *
 * The Hall sets will almost always be instantiated as part of the Lie basis,
 * which is usually a static member of a vector. For this reason, it makes it
 * tricky to include these data as part of the Hall set object itself. Instead
 * we make use of a container class which can be instantiated as a static variable
 * in a function, and thus be instantiated only once per width, and the main Hall
 * set class can get store a reference to this static object. This model also
 * greatly simplifies the growing policy, since we can build in the grow_up call
 * into the function that generates the instance.
 *
 * In terms of content, this class contains seven individual members: a vector
 * of letters; a vector of parents (key to parents); a reverse map from parents
 * to keys; the current degree held; a map from letters to keys; a vector of
 * degree ranges; a vector of sizes (size of the Hall set up to degree).
 * These are all public members.
 *
 * @tparam Width The alphabet size of the Hall set.
 */
template<DEG Width>
class hall_set_content {
public:
    using degree_type = DEG;
    using letter_type = LET;
    using key_type = letter_type;
    using size_type = DIMN;
    using parent_type = std::pair<key_type, key_type>;

    using degree_range_type = std::pair<size_type, size_type>;
    using reverse_map_type = std::map<parent_type, key_type>;
    using data_type = std::vector<parent_type>;
    using letter_vec_type = std::vector<letter_type>;
    using l2k_map_type = std::vector<key_type>;
    using degree_ranges_map_type = std::vector<degree_range_type>;
    using size_map_type = std::vector<size_type>;

    letter_vec_type letters;
    data_type data;
    reverse_map_type reverse_map;
    degree_type current_degree;
    l2k_map_type l2k;
    degree_ranges_map_type degree_ranges;
    std::vector<size_type> sizes;

    static constexpr degree_type n_letters = Width;

private:

    hall_set_content() : current_degree{0}
    {
        data.reserve(1+n_letters);
        sizes.reserve(2);
        l2k.reserve(n_letters);
        data.push_back(parent_type(0, 0));
        degree_ranges.reserve(2);
        degree_ranges.push_back(degree_range_type(0, 1));
        sizes.push_back(0);

        for (letter_type l = 1; l<=n_letters; ++l) {
            parent_type parents(key_type(0), key_type(l));
            letters.push_back(l);
            data.push_back(parents);
            reverse_map.insert(std::pair<parent_type, key_type>(parents, data.size()-1));
            l2k.push_back(l);
        }

        degree_range_type range;
        range.first = degree_ranges[current_degree].second;
        range.second = data.size();
        degree_ranges.push_back(range);
        sizes.push_back(n_letters);
        ++current_degree;
    }

    void grow_up(degree_type degree)
    {
        for (degree_type d = current_degree+1; d<=degree; ++d) {
            for (degree_type e = 1; 2*e<=d; ++e) {
                letter_type i_lower, i_upper, j_lower, j_upper;
                i_lower = degree_ranges[e].first;
                i_upper = degree_ranges[e].second;
                j_lower = degree_ranges[d-e].first;
                j_upper = degree_ranges[d-e].second;

                for (letter_type i = i_lower; i<i_upper; ++i) {
                    for (letter_type j = std::max(j_lower, i+1); j<j_upper; ++j) {
                        if (data[j].first<=i) {
                            parent_type parents(i, j);
                            data.push_back(parents);
                            reverse_map[parents] = data.size()-1;
                        }
                    }
                }
            }

            degree_range_type range;
            range.first = degree_ranges[current_degree].second;
            range.second = data.size();
            degree_ranges.push_back(range);

            ++current_degree;
        }
    }

public:

    static hall_set_content& instance(degree_type degree = 0) noexcept
    {
        static std::mutex lock;
        std::lock_guard<std::mutex> access(lock);

        static hall_set_content content;
        if (degree>=2) {
            content.grow_up(degree);
        }

        return content;
    }

};


template <DEG Width>
class hall_set
{

public:
    using degree_type               = DEG;
    using letter_type               = LET;
    using key_type                  = letter_type;
    using size_type                 = DIMN;
    using parent_type               = std::pair<key_type, key_type>;

private:
    using hall_set_content_type     = hall_set_content<Width>;
    using degree_range_type         = std::pair<size_type, size_type>;
    using reverse_map_type          = std::map<parent_type, key_type>;
    using data_type                 = std::vector<parent_type>;
    using letter_vec_type           = std::vector<letter_type>;
    using l2k_map_type              = std::vector<key_type>;
    using degree_ranges_map_type    = std::vector<degree_range_type>;
    using size_map_type             = std::vector<size_type>;

    hall_set_content_type& content;

public:

    static constexpr degree_type n_letters = Width;

private:

    hall_set()
        : content(hall_set_content_type::instance()),
          key2string(*this, letter_to_string, key2string_binop)
    {
    }

public:

    explicit hall_set(degree_type degree)
        :content(hall_set_content_type::instance(degree)),
          key2string(*this, letter_to_string, key2string_binop)
    {
    }


public:

    const letter_vec_type& letters() const noexcept
    {
        return content.letters;
    }

    const degree_ranges_map_type& hall_set_degree_ranges() const noexcept
    {
        return content.degree_ranges;
    }

    const l2k_map_type& l2k() const noexcept
    {
        return content.l2k;
    }

private:

    using signed_size_t = std::make_signed<size_type>::type;

    static constexpr signed_size_t
    level_term(signed_size_t degree, signed_size_t divisor)
    {
        return mobius(divisor)*power(static_cast<signed_size_t>(Width), degree/divisor) / degree;
    }

    static constexpr signed_size_t
    level_size_impl(signed_size_t degree, signed_size_t divisor)
    {
        return ((degree % divisor == 0) ? level_term(degree, divisor) : 0)
                + ((divisor < degree) ? level_size_impl(degree, divisor+1) : 0);
    }

    static constexpr size_type level_size(degree_type degree) noexcept
    {
        return static_cast<size_type>(level_size_impl(degree, 1));
    }

public:

    static constexpr size_type start_of_degree(degree_type degree) noexcept
    {
        return (degree == 0 || degree == 1) ? 0 : (level_size(degree-1) + start_of_degree(degree-1));
    }

    key_type begin() const noexcept
    {
        return key_type(1);
    }

    key_type end() const noexcept
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


    key_type keyofletter(letter_type let) const
    {
        assert(letter(let));
        return content.l2k[let-1];
    }

    key_type lparent(const key_type& key) const noexcept
    {
        return content.data[key].first;
    }

    key_type rparent(const key_type& key) const noexcept
    {
        return content.data[key].second;
    }

    bool letter(const key_type& key) const noexcept
    {
        return ((key > 0) && (key <= Width));
    }

    letter_type getletter(const key_type& key) const
    {
        return content.letters[key - 1];
    }

    size_type size() const noexcept
    {
        return content.data.size() - 1;
    }

    size_type key_to_index(const key_type& key) const noexcept
    {
        return static_cast<size_type>(key - 1);
    }

    key_type index_to_key(size_type index) const noexcept
    {
        return static_cast<key_type>(index + 1);
    }

    const parent_type& operator[](const key_type& key) const noexcept
    {
        return content.data[key];
    }

    const key_type& operator[](const parent_type& parents) const
    {
        typename reverse_map_type::const_iterator it;
        if ((it = content.reverse_map.find(parents)) != content.reverse_map.end()) {
            return it->second;
        } else {
            throw std::invalid_argument("parent object does no correspond to hall basis element");
        }
    }

    typename reverse_map_type::const_iterator find(const parent_type& parents) const noexcept
    {
        return content.reverse_map.find(parents);
    }

    typename reverse_map_type::const_iterator reverse_map_end() const noexcept
    {
        return content.reverse_map.end();
    }

    typename data_type::const_iterator parents_begin() const noexcept
    {
        return content.data.cbegin();
    }

    typename data_type::const_iterator parents_end() const noexcept
    {
        return content.data.cend();
    }

public:


    template<typename Function, typename BinOp, typename Tag>
    class extended_function {
    public:
        typedef Function function_type;
        typedef BinOp binary_operation_type;
        typedef Tag tag_type;
        typedef decltype(std::declval<function_type>()(std::declval<key_type>())) output_type;

        using table_t = std::map<key_type, output_type>;

        extended_function(const hall_set& hs)
                : m_hall_set(hs), m_tag(), m_fn(), m_op() { }

        extended_function(const hall_set& hs, Function fn, BinOp op)
                : m_hall_set(hs), m_tag(), m_fn(fn), m_op(op) { }

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

        template<typename Predicate>
        output_type eval(const key_type& k, lazy_cache_tag<Predicate> tag) const
        {
            static boost::recursive_mutex table_lock;
            static table_t table;

            if (tag.predicate(k)) {
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

        output_type eval(const key_type& k, lazy_cache_tag<void>) const
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

        template<typename Predicate>
        table_t fill_table(Predicate predicate) const
        {
            table_t result;

            for (key_type k=m_hall_set.begin(); k!=m_hall_set.end() && predicate(k); k = m_hall_set.next_key(k)) {
                if (m_hall_set.letter(k)) {
                    result[k] = m_fn(m_hall_set.getletter(k));
                } else {
                    result[k] = m_op(result[m_hall_set.lparent(k)], result[m_hall_set.rparent(k)]);
                }
            }

            return result;
        }

        template<typename Predicate>
        output_type eval(const key_type& k, lookup_table_tag<Predicate>) const
        {
            static table_t table = fill_table(m_tag.predicate);

            if (m_tag.predicate(k)) {
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
        const hall_set& m_hall_set;
        function_type m_fn;
        binary_operation_type m_op;
        tag_type m_tag;
    };


    using key2string_type = extended_function<
            decltype(&hall_set::letter_to_string),
            decltype(&hall_set::key2string_binop),
            lazy_cache_tag<void>
    >;

    const key2string_type key2string;

};






}




#endif //LIBALGEBRA_HALL_SET_H
