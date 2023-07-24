//
// Created by sam on 23/09/2021.
//

#ifndef LIBALGEBRA_HALL_SET_H
#define LIBALGEBRA_HALL_SET_H

#include <cassert>
#include <map>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#endif

#include "implementation_types.h"
#include "detail/integer_maths.h"
#include "detail/caching_tags.h"
#include "detail/function_extension_cache_base.h"

namespace alg {

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
class hall_set_content
{
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
    hall_set_content()
        : current_degree{0}
    {
        data.reserve(1 + n_letters);
        letters.reserve(n_letters);
        sizes.reserve(2);
        l2k.reserve(n_letters);
        data.push_back(parent_type(0, 0));
        degree_ranges.reserve(2);
        degree_ranges.push_back(degree_range_type(0, 1));
        sizes.push_back(0);

        for (letter_type l = 1; l <= n_letters; ++l) {
            parent_type parents(key_type(0), key_type(l));
            letters.push_back(l);
            data.push_back(parents);
            reverse_map.insert(std::pair<parent_type, key_type>(parents, data.size() - 1));
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
        for (degree_type d = current_degree + 1; d <= degree; ++d) {
            for (degree_type e = 1; 2 * e <= d; ++e) {
                letter_type i_lower, i_upper, j_lower, j_upper;
                i_lower = degree_ranges[e].first;
                i_upper = degree_ranges[e].second;
                j_lower = degree_ranges[d - e].first;
                j_upper = degree_ranges[d - e].second;

                for (letter_type i = i_lower; i < i_upper; ++i) {
                    for (letter_type j = std::max(j_lower, i + 1); j < j_upper; ++j) {
                        if (data[j].first <= i) {
                            parent_type parents(i, j);
                            data.push_back(parents);
                            reverse_map[parents] = data.size() - 1;
                        }
                    }
                }
            }

            degree_range_type range;
            range.first = degree_ranges[current_degree].second;
            range.second = data.size();
            degree_ranges.push_back(range);
            sizes.push_back(data.size());

            ++current_degree;
        }
    }

public:
    static hall_set_content& instance(degree_type degree = 0)
    {
        static std::mutex lock;
        std::lock_guard<std::mutex> access(lock);

        static hall_set_content content;
        if (degree >= 2) {
            content.grow_up(degree);
        }

        return content;
    }
};

template <DEG Width>
constexpr typename hall_set_content<Width>::degree_type hall_set_content<Width>::n_letters;


/**
 * @brief A Hall set
 *
 * A basis is a finite total ordered set of keys, its cardinal is size() and
 * its minimal element is begin(). The successor key of a given key is given
 * by nextkey(). The successor of the maximal key is end() and does not belong
 * to the basis. The position of a given key in the total order of the basis
 * is given by keypos(), and equals 1 for begin(). To each letter corresponds
 * a key.\n\n
 *
 * This class is an ancestor of the lie_basis class, used to implement the lie
 * class (Free Lie Algebra) as a particular instance of an algebra class
 * (Associative Algebras).\n\n
 *
 * This class stores a Philip Hall basis associated to a finite number of
 * letters. A key is the implementation of a Lie element of this basis. A
 * letter is a particular Lie element (or basis element, or key). Each key k
 * which does not correspond to a letter has two parents lp and rp and we have
 * k = [lp,rp] where [.,.] is the Lie product. A letter, viewed as a key, has
 * no parents. More precisely, its parents are invalid keys.\n\n
 *
 * The basis elements are recursively computed and are enumerated with keys.
 * The set of valid keys is essentially an interval of natural integers.\n\n
 *
 * One can find below a brief Mathematical description of Philip Hall bases
 * for the free Lie Algebra. Cf. Reutenauer's book for example, ISBN 0 19
 * 853679 8.\n\n
 *
 * Let K be a field with characteristic non equals to 2. In
 * newgenesis-libalgebra, this field K corresponds to the type SCA defined in
 * libalgebra.h.\n\n
 *
 * Let M be a finite alphabet {a_1,...,a_n}. We denote by M* the monoid which
 * consists in words of letters in M. The product in M* is the concatenation
 * and the neutral element is the empty word.\n\n
 *
 * We consider the free algebra A over (K,M). An element of A is a linear
 * combination of elements of M*, with coefficients in K. An element of A is
 * an instance of class free_tensor, which affects to each element of M* a
 * coefficient in K. The element of M* are indexed by tensor_key, which
 * essentially stores the corresponding word as a std::string.\n\n
 *
 * We consider also the associated free Lie algebra L, the smallest subalgebra
 * of A which contains M and is stable by the Lie product [X,Y] = XY-YX. An
 * element of L is an instance of class lie. The key used are of type
 * lie_key, which are actually indexes in a basis of type lie_basis.\n\n
 *
 * The degree of a word w in M is its length. The degree of an element of the
 * algebra A is the maximum degree of words with non-zero coefficients. The
 * degree of [X,Y] is the sum of the degrees of X and Y if X and Y are
 * different, and 0 if X = Y.\n\n
 *
 * Actually, the free Lie algebra L is a graded algebra, with respect to the
 * degree (or weight) of Lie products. Philip Hall invented an algorithm for
 * computing a basis of the free Lie algebra L. A Hall basis H is a union of
 * subsets H_1,...H_i,... of L. By definition, H_1 = M = {a_1,...,a_n} and the
 * elements of H_i are of degree i. The set H is totally ordered and more
 * over, H_1 \< H_2 \< ... The Hall basis H can be constructed recursively from
 * H_1. This can be done by constructing an array HALLARRAY of elements of the
 * form {left, degree, right}. The left and right corresponds to indexes in
 * the array for constructing the element by the Lie product, and degree
 * corresponds to the degree of this element, which is then the sum of the
 * degrees of the two elements pointed by left and right. The order of the
 * elements of the array is in one to one correspondence with the order of H.
 * The subset H_i is exactly the elements of the form {left, degree, right}
 * with degree = i.\n\n
 *
 * Starting from H1 = {{0, 1, 1},...,{0, 1, n}} which corresponds to the n
 * letters, Hi+1 is constructed from H_1, ..., H_i by examining all elements
 * of the form {l, i + 1, r} where l < r and l and r are in the union of
 * H_1,...,H_i. Such an element is added to the set Hi+1 if and only if the
 * right parent of r is \<= l.
 *
 * @tparam Width Size of alphabet
 */
template<DEG Width>
class hall_set
{

public:
    using degree_type = DEG;
    using letter_type = LET;
    using key_type = letter_type;
    using size_type = DIMN;
    using parent_type = std::pair<key_type, key_type>;

protected:
    using hall_set_content_type = hall_set_content<Width>;
    using degree_range_type = std::pair<size_type, size_type>;
    using reverse_map_type = std::map<parent_type, key_type>;
    using data_type = std::vector<parent_type>;
    using letter_vec_type = std::vector<letter_type>;
    using l2k_map_type = std::vector<key_type>;
    using degree_ranges_map_type = std::vector<degree_range_type>;
    using size_map_type = std::vector<size_type>;

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
        : content(hall_set_content_type::instance(degree)),
          key2string(*this, letter_to_string, key2string_binop)
    {
    }

public:
    /// Get a const reference to vector of letters
    const letter_vec_type& letters() const noexcept
    {
        return content.letters;
    }

    /// Get a const reference to vector of degree ranges
    const degree_ranges_map_type& hall_set_degree_ranges() const noexcept
    {
        return content.degree_ranges;
    }

    const degree_type& current_degree() const noexcept
    {
        return content.current_degree;
    }

    /// return a const reference to letter to key map
    const l2k_map_type& l2k() const noexcept
    {
        return content.l2k;
    }

private:
    using signed_size_t = std::make_signed<size_type>::type;

    static constexpr signed_size_t
    level_term(signed_size_t degree, signed_size_t divisor)
    {
        return mobius(divisor) * power(static_cast<signed_size_t>(Width), degree / divisor) / degree;
    }

    static constexpr signed_size_t
    level_size_impl(signed_size_t degree, signed_size_t divisor)
    {
        return ((degree % divisor == 0) ? level_term(degree, divisor) : 0)
                + ((divisor < degree) ? level_size_impl(degree, divisor + 1) : 0);
    }

    static constexpr size_type level_size(degree_type degree) noexcept
    {
        return static_cast<size_type>(level_size_impl(degree, 1));
    }

public:
    /// Get the index in the basis at which the elements of given degree start
    static constexpr size_type start_of_degree(degree_type degree) noexcept
    {
        return (degree == 0 || degree == 1) ? 0 : (level_size(degree - 1) + start_of_degree(degree - 1));
    }

    /// Get the first element of the Hall set
    key_type begin() const noexcept
    {
        return key_type(1);
    }

    /// Get the "end" of the hall set - this is the invalid key 0
    key_type end() const noexcept
    {
        return key_type(0);
    }

    /// Get the key that follows in the basis order
    key_type next_key(const key_type& key) const noexcept
    {
        if (key < size()) {
            return (key + 1);
        }
        else {
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
    /// Get the key that corresponds to a letter
    key_type keyofletter(letter_type let) const
    {
        assert(letter(let));
        return content.l2k[let - 1];
    }

    /// Get the left parent of a key
    key_type lparent(const key_type& key) const noexcept
    {
        return content.data[key].first;
    }

    /// Get the right parent of a key
    key_type rparent(const key_type& key) const noexcept
    {
        return content.data[key].second;
    }

    /// Test if a key corresponds to a letter
    bool letter(const key_type& key) const noexcept
    {
        return ((key > 0) && (key <= Width));
    }

    /// Get the letter that corresponds to a key.
    letter_type getletter(const key_type& key) const
    {
        assert(letter(key));
        return content.letters[key - 1];
    }

    /// Get the current size of the Hall set
    size_type size() const noexcept
    {
        return content.data.size() - 1;
    }

    /// Get the index at which a key appears in the basis order
    size_type key_to_index(const key_type& key) const noexcept
    {
        return static_cast<size_type>(key - 1);
    }

    /// Get the key that appears at a given index
    key_type index_to_key(size_type index) const noexcept
    {
        return static_cast<key_type>(index + 1);
    }

    /// Retrieve the parents of a key
    const parent_type& operator[](const key_type& key) const noexcept
    {
        return content.data[key];
    }

    /// Retrieve the key that is formed from the bracket of parents - throws an exception if invalid
    const key_type& operator[](const parent_type& parents) const
    {
        typename reverse_map_type::const_iterator it;
        if ((it = content.reverse_map.find(parents)) != content.reverse_map.end()) {
            return it->second;
        }
        else {
            throw std::invalid_argument("parent object does no correspond to hall basis element");
        }
    }

    /// Retrieve the key that is formed from the bracket of parents
    typename reverse_map_type::const_iterator find(const parent_type& parents) const noexcept
    {
        return content.reverse_map.find(parents);
    }

    /// Get the iterator end of the reverse map
    typename reverse_map_type::const_iterator reverse_map_end() const noexcept
    {
        return content.reverse_map.end();
    }

    /// Get the iterator at the start of the Hall set data
    typename data_type::const_iterator parents_begin() const noexcept
    {
        return content.data.cbegin();
    }

    /// Get the pointer to the end of the Hall set data
    typename data_type::const_iterator parents_end() const noexcept
    {
        return content.data.cend();
    }

public:
    /**
     * @brief A wrapper around a function from letters into some set with a binary operation
     *
     * This is a wrapper around a function-like type and a binary operation type that is used
     * to extend a function to the Hall set by recursively computing the result of
     * f([k1, k2]) = f(k1)*f(k2), where f is the extended function and * is the binary operation.
     * This recursion eventually terminates on letters, where the function is initially defined,
     * so this process works and results in the desired extension.
     *
     * @tparam Function Function type
     * @tparam BinOp Binary operation type
     * @tparam Tag Cache type indicator tag
     * @see hall_basis::extend_function
     */
    template<typename Function, typename BinOp, typename Tag>
    class extended_function
    {
    public:
        typedef Function function_type;
        typedef BinOp binary_operation_type;
        typedef Tag tag_type;
        typedef decltype(std::declval<function_type>()(std::declval<key_type>())) output_type;

        using table_t = std::map<key_type, output_type>;

        // Default constructors
        extended_function() = default;
        extended_function(const extended_function&) = default;
        extended_function(extended_function&&) noexcept = default;

        explicit extended_function(const hall_set& hs)
            : m_hall_set(hs), m_fn(), m_op(), m_tag()
        {}

        template<typename HallSet>
        explicit extended_function(const HallSet& hs) : m_hall_set(hs), m_fn(), m_op(), m_tag()
        {}

        template<typename HallSet>
        extended_function(const HallSet& hs, Function fn) : m_hall_set(hs), m_fn(fn), m_op(), m_tag()
        {}

        template<typename HallSet>
        extended_function(const HallSet& hs, BinOp bin_op) : m_hall_set(hs), m_fn(), m_op(bin_op), m_tag()
        {}

        template<typename HallSet>
        extended_function(const HallSet& hs, Function fn, BinOp bin_op) : m_hall_set(hs), m_fn(fn), m_op(bin_op), m_tag()
        {}

        extended_function(const hall_set& hs, Function fn, BinOp op)
            : m_hall_set(hs), m_fn(fn), m_op(op), m_tag()
        {}

    private:
        output_type eval_impl(const key_type& k) const
        {
            if (m_hall_set.letter(k)) {
                return m_fn(m_hall_set.getletter(k));
            }
            else {
                return m_op(
                        operator()(m_hall_set.lparent(k)),
                        operator()(m_hall_set.rparent(k)));
            }
        }

        output_type eval(const key_type& k, no_caching_tag) const
        {
            return eval_impl(k);
        }

        template<typename Predicate>
        output_type eval(const key_type& k, lazy_cache_tag<Predicate> tag) const
        {
            static std::recursive_mutex table_lock;
            static table_t table;

            if (tag.predicate(k)) {
                std::lock_guard<std::recursive_mutex> access(table_lock);

                typename table_t::iterator it = table.find(k);
                if (it != table.end()) {
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
            static std::recursive_mutex table_lock;
            static table_t table;

            std::lock_guard<std::recursive_mutex> access(table_lock);

            typename table_t::iterator it = table.find(k);
            if (it != table.end()) {
                return it->second;
            }

            return table[k] = eval_impl(k);
        }

        template<typename Predicate>
        table_t fill_table(Predicate predicate) const
        {
            table_t result;

            for (key_type k = m_hall_set.begin(); k != m_hall_set.end() && predicate(k); k = m_hall_set.next_key(k)) {
                if (m_hall_set.letter(k)) {
                    result[k] = m_fn(m_hall_set.getletter(k));
                }
                else {
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

        output_type eval(const key_type& k, lazy_cache_on_object_tag) const
        {
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
            lazy_cache_tag<void>>;

    /// Write out a key as a (nested) bracket of its parents
    const key2string_type key2string;

#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
protected:
    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive& ar, const unsigned int /*version*/)
    {
        ar& content.letters;
        ar& content.data;
        ar& content.reverse_map;
        ar& content.current_degree;
        ar& content.l2k;
        ar& content.degree_ranges;
        ar& content.sizes;
    }
#endif
};

}// namespace alg

#endif//LIBALGEBRA_HALL_SET_H
