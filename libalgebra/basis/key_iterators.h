//
// Created by sam on 16/11/2021.
//

#ifndef LIBALGEBRA_KEY_ITERATORS_H
#define LIBALGEBRA_KEY_ITERATORS_H

#include <iterator>

namespace alg {
namespace basis {

template<typename Basis>
class basis_iterable;

template<typename Basis>
class key_iterator
{
    using parent_type = basis_iterable<Basis>;

public:
    using difference_type = std::ptrdiff_t;
    using value_type = typename Basis::KEY;
    using pointer = const value_type*;
    using reference = const value_type&;
    using iterator_category = std::forward_iterator_tag;

    key_iterator() : parent(nullptr), current()
    {}

    explicit key_iterator(parent_type* p) : parent(p), current(p->basis.begin())
    {}

    key_iterator(parent_type* p, value_type c) : parent(p), current(c)
    {}

    key_iterator(const key_iterator& other) : parent(other.parent), current(other.current)
    {}

    key_iterator(key_iterator&& other) noexcept
    {
        parent = other.parent;
        current = std::move(other.current);
    }

    key_iterator& operator=(const key_iterator& other)
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

    reference operator*()
    {
        return current;
    }

    pointer operator->()
    {
        return &current;
    }

public:
    bool operator==(const key_iterator& other) const
    {
        return (parent == other.parent && current == other.current);
    }

    bool operator!=(const key_iterator& other) const
    {
        return (parent != other.parent || current != other.current);
    }

private:
    parent_type* parent;
    value_type current;
};

template<typename Basis>
class basis_iterable
{
    friend class key_iterator<Basis>;

public:
    using iterator = key_iterator<Basis>;
    using value_type = typename Basis::KEY;

    explicit basis_iterable(const Basis& b) : basis(b), start_key(basis.begin()), end_key(basis.end())
    {}

    explicit basis_iterable(const Basis& b, const value_type& start) : basis(b), start_key(start), end_key(basis.end())
    {}

    explicit basis_iterable(const Basis& b, const value_type& start, const value_type& end)
        : basis(b), start_key(start), end_key(basis.nextkey(end))
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
