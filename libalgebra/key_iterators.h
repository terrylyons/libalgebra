//
// Created by sam on 16/11/2021.
//

#ifndef LIBALGEBRA_KEY_ITERATORS_H
#define LIBALGEBRA_KEY_ITERATORS_H

#include <iterator>

namespace alg {
namespace basis {

template<typename Basis>
class key_range;

/**
 * @brief Iterator of keys of a basis.
 *
 * A key_iterator of a basis produces the keys associated with the basis in order.
 * This acts in tandem with a parent object that controls the start and end of
 * iteration. Two iterators with different parents will always be non-equal.
 * This is a forward input iterator - meaning it cannot be used to modify the
 * keys of a basis. Indeed, the keys in a basis must be invariant.
 *
 * Under the hood, the key_iterator uses the basis function nextkey to advance
 * the iterator, and holds the current key internally.
 *
 * @tparam Basis Basis type
 */
template<typename Basis>
class key_iterator
{
    using parent_type = key_range<Basis>;

public:
    using difference_type = std::ptrdiff_t;
    using value_type = typename Basis::KEY;
    using pointer = const value_type*;
    using reference = const value_type&;
    using iterator_category = std::forward_iterator_tag;

    key_iterator() noexcept : parent(nullptr), current()
    {}

    explicit key_iterator(parent_type* p) noexcept : parent(p), current(p->basis.begin())
    {}

    key_iterator(parent_type* p, value_type c) noexcept : parent(p), current(c)
    {}

    key_iterator(const key_iterator& other) noexcept : parent(other.parent), current(other.current)
    {}

    key_iterator(key_iterator&& other) noexcept
    {
        parent = other.parent;
        current = std::move(other.current);
    }

    key_iterator& operator=(const key_iterator& other) noexcept
    {
        parent = other.parent;
        current = other.current;
        return *this;
    }

    key_iterator& operator=(key_iterator&& other) noexcept
    {
        parent = other.parent;
        current = std::move(other.current);
        return *this;
    }

public:
    key_iterator& operator++()
    {
        current = parent->basis.nextkey(current);
        return *this;
    }

    key_iterator operator++(int)
    {
        key_iterator prev(parent, current);
        ++(*this);
        return prev;
    }

    reference operator*() noexcept
    {
        return current;
    }

    pointer operator->() noexcept
    {
        return &current;
    }

public:
    bool operator==(const key_iterator& other) const noexcept
    {
        return (parent == other.parent && current == other.current);
    }

    bool operator!=(const key_iterator& other) const noexcept
    {
        return (parent != other.parent || current != other.current);
    }

private:
    parent_type* parent;
    value_type current;
};

/**
 * @brief Iterable class representing a range of keys associated with a basis.
 *
 * This class allows us to iterate over the keys associated with a basis from a given range.
 * The class holds a reference to the basis, and so must be provided with this reference
 * upon construction. By default, when only this basis reference is provided, the range will
 * contain all keys associated with the basis. Alternatively, one can provide a start key
 * or both start key and end key. In the first case, the range will iterate over all keys
 * from the start key provided to the end, and in the second case the range will iterate
 * from start key (inclusive) to end key (exclusive).
 *
 * The iterator type of the range is a key_iterator<Basis>.
 *
 * @tparam Basis Basis type from which keys should be taken
 */
template<typename Basis>
class key_range
{
    friend class key_iterator<Basis>;

public:
    using iterator = key_iterator<Basis>;
    using value_type = typename Basis::KEY;

    explicit key_range(const Basis& b) : basis(b), start_key(basis.begin()), end_key(basis.end())
    {}

    explicit key_range(const Basis& b, const value_type& start) : basis(b), start_key(start), end_key(basis.end())
    {}

    explicit key_range(const Basis& b, const value_type& start, const value_type& end)
        : basis(b), start_key(start), end_key(end)
    {}

    iterator begin() noexcept
    {
        return iterator(this, start_key);
    }

    iterator end() noexcept
    {
        return iterator(this, end_key);
    }

private:
    const Basis& basis;
    value_type start_key;
    value_type end_key;
};

}// namespace basis
}// namespace alg

#endif//LIBALGEBRA_KEY_ITERATORS_H
