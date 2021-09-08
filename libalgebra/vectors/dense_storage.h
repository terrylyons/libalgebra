//
// Created by sam on 04/08/2021.
//

#ifndef LIBALGEBRA_DENSE_STORAGE_H
#define LIBALGEBRA_DENSE_STORAGE_H

#include <memory>
#include <cassert>
#include <iostream>

#include "libalgebra/implimentation_types.h"

namespace alg {
namespace vectors {


namespace dtl {

template <typename S, typename Alloc>
struct dense_storage_base
{
    using allocator_type = Alloc;
    using alloc_traits   = std::allocator_traits<allocator_type>;
    using value_type     = S;
    using pointer        = S*;
    using const_pointer  = S const*;
    using size_type      = typename allocator_type::size_type;

    enum vec_type
    {
        owned, borrowed_mut, borrowed
    };

    vec_type m_type;
    pointer m_data;
    size_type m_size;
    allocator_type m_alloc;

    explicit dense_storage_base(size_type sz=0)
        : m_alloc{},
        m_data{(sz > 0) ? alloc_traits::allocate(m_alloc, sz) : nullptr},
        m_size{sz}, m_type{owned}
    {
    }

    ~dense_storage_base()
    {
        if (m_type == owned) {
            alloc_traits::deallocate(m_alloc, m_data, m_size);
        }
    }

    dense_storage_base(pointer begin, pointer end)
        : m_alloc{}, m_data{begin}, m_size{static_cast<size_type>(end - begin)}, m_type{borrowed_mut}
    {}

    dense_storage_base(const_pointer begin, const_pointer end)
        : m_alloc{},
          m_data{const_cast<pointer>(begin)},
          m_size{static_cast<size_type>(end - begin)},
          m_type{borrowed}
    {}

    dense_storage_base(dense_storage_base const&) = delete;
    dense_storage_base& operator=(dense_storage_base const&) = delete;

    dense_storage_base(dense_storage_base&& other) noexcept
        : m_alloc{other.m_alloc}, m_data{other.m_data}, m_size{other.m_size}, m_type{other.m_type}
    {
        other.m_data = nullptr;
        other.m_size = 0;
        other.m_type = owned;
    }

    dense_storage_base& operator=(dense_storage_base&& other) noexcept
    {
        std::swap(m_alloc, other.m_alloc);
        std::swap(m_data, other.m_data);
        std::swap(m_type, other.m_type);
        std::swap(m_size, other.m_size);
        return *this;
    }


    bool is_owned() const
    {
        return m_type == owned;
    }

    bool is_borrowed() const
    {
        return m_type == borrowed;
    }

    bool is_borrowed_mut() const
    {
        return m_type == borrowed_mut;
    }

};


}


/**
 * @brief Storage class for dense vectors.
 *
 * This class will handle the storage for the dense vector type. For the time being, this will simply
 * hold a pointer to the data, the size of the buffer and the type of the vector which can be either
 * `owned`, `borrowed_mut`, `borrowed`. This storage type will implement a "copy on resize" for `borrowed`
 * and `borrowed_mut` types, and "copy on modify" for "borrowed" types.
 *
 * For the time being, this is going to use `new` and `delete` to manage it's memory, which I know is bad.
 * However, I want this to be able to allocate without assigning since we can potentially waste a lot of
 * filling a vector only to then replace all the entries by assignment.
 *
 *
 * @tparam S
 */
template <typename S, typename Alloc = std::allocator<S> >
class dense_storage
{
    using base_type = dtl::dense_storage_base<S, Alloc>;
public:
    using allocator_type = typename base_type::allocator_type;
    using alloc_traits   = typename base_type::alloc_traits;

    using size_type = typename base_type::size_type;

    using is_pod_t = std::is_pod<S>;

    using value_type = S;
    using reference = S &;
    using const_reference = S const &;
    using pointer = S *;
    using const_pointer = S const *;

    using iterator = S *;
    using const_iterator = S const *;

    using vec_type = typename base_type::vec_type;

private:

    dtl::dense_storage_base<value_type, allocator_type> m_base;

    static void destroy_range(pointer, pointer, std::true_type)
    {
    }

    static void destroy_range(pointer start, pointer end, std::false_type)
    {
        while (end != start) {
            --end;
            end->~S();
        }
    }

    static void destroy_range(pointer start, pointer end)
    {
        std::is_trivially_default_constructible<S> tag;
        destroy_range(start, end, tag);
    }

    void fill_range_default_construct(pointer start, pointer end)
    {
        std::uninitialized_fill(start, end, value_type());
    }


public:


    explicit dense_storage(size_type size = 0) : m_base{size}
    {
        fill_range_default_construct(m_base.m_data, m_base.m_data + m_base.m_size);
    }

    dense_storage(dense_storage const &other) : m_base{other.size()}
    {
        std::uninitialized_copy(other.begin(), other.end(), m_base.m_data);
    }

    dense_storage(pointer ptr, size_type size) : m_base{ptr, ptr + size}
    {
    }

    dense_storage(const_pointer ptr, size_type size) : m_base{ptr, ptr + size}
    {
    }

    dense_storage(pointer begin, pointer end) : m_base{begin, end}
    {
    }

    dense_storage(const_pointer begin, const_pointer end) : m_base{begin, end}
    {
    }

    dense_storage(size_type offset, const_pointer start, const_pointer end)
        : m_base{offset + static_cast<size_type>(end - start)}
    {
        fill_range_default_construct(m_base.m_data, m_base.m_data + offset);
        std::uninitialized_copy(start, end, m_base.m_data+offset);
    }

    dense_storage(size_type offset, pointer start, pointer end)
        : m_base{offset + static_cast<size_type>(end - start)}
    {
        fill_range_default_construct(m_base.m_data, m_base.m_data + offset);
        std::uninitialized_copy(std::make_move_iterator(start), std::make_move_iterator(end), m_base.m_data + offset);
    }

    ~dense_storage()
    {
        if (m_base.is_owned()) {
            destroy_range(m_base.m_data, m_base.m_data+m_base.m_size);
        }
    }

    dense_storage &operator=(dense_storage const &other)
    {
        dense_storage tmp(other);
        this->swap(tmp);
        return *this;
    }

    dense_storage &operator=(dense_storage &&other)
    {
        if (m_base.is_owned()) {
            destroy_range(m_base.m_data, m_base.m_data + m_base.m_size);
        }
        m_base = std::move(other.m_base);
        return *this;
    }


    static dense_storage make_owned(const_pointer ptr, size_type sz)
    {
        dense_storage result(ptr, sz);
        result.to_owned();
        return result;
    }


private:

    void maybe_fill(pointer, pointer, const_reference, std::true_type)
    {}

    void maybe_fill(pointer range_begin, pointer range_end, const_reference val, std::false_type)
    {
        std::uninitialized_fill(range_begin, range_end, val);
    }

    void maybe_fill(pointer range_begin, pointer range_end, const_reference val)
    {
        std::is_trivially_default_constructible<value_type> tag;
        maybe_fill(range_begin, range_end, val, tag);
    }

    void to_owned(size_type alloc_size)
    {
        assert(!m_base.is_owned());
        base_type new_base(alloc_size);
        size_type mid = std::min(alloc_size, size());

        std::uninitialized_copy(m_base.m_data, m_base.m_data+mid, new_base.m_data);
        if (alloc_size > mid) {
            maybe_fill(new_base.m_data+mid, new_base.m_data+new_base.m_size, value_type());
        }

        m_base = std::move(new_base);
    }

    void to_owned(size_type alloc_size, const_reference val)
    {
        assert(!m_base.is_owned());
        base_type new_base(alloc_size);
        size_type mid = std::min(alloc_size, size());

        std::uninitialized_copy(m_base.m_data, m_base.m_data+m_base.m_size, new_base.m_data);

        if (mid < alloc_size) {
            // new larger than old
            std::uninitialized_fill(new_base.m_data+m_base.m_size, new_base.m_data+alloc_size, val);
        }

        m_base = std::move(new_base);
    }

    void resize_owned(size_type sz)
    {
        assert(m_base.is_owned());
        assert(size() != sz);

        size_type mid = std::min(sz, size());

        base_type new_base(sz);
        std::uninitialized_copy(
                std::make_move_iterator(m_base.m_data),
                std::make_move_iterator(m_base.m_data+mid),
                new_base.m_data
                );

        if (mid >= sz) {
            destroy_range(m_base.m_data + mid, m_base.m_data+m_base.m_size);
        } else {
            maybe_fill(new_base.m_data+mid, new_base.m_data+new_base.m_size, value_type());
        }

        m_base = std::move(new_base);
    }

    void resize_owned(size_type sz, const_reference val)
    {
        assert(m_base.is_owned());
        assert(size() != sz);

        size_type mid = std::min(sz, size());

        base_type new_base(sz);
        std::uninitialized_copy(
                std::make_move_iterator(m_base.m_data),
                std::make_move_iterator(m_base.m_data+mid),
                new_base.m_data
                );

        if (mid < size()) {
            destroy_range(m_base.m_data + mid, m_base.m_data+m_base.m_size);
        } else {
            std::uninitialized_fill(new_base.m_data + mid, new_base.m_data + sz, val);
        }

        m_base = std::move(new_base);
    }

    void to_owned()
    {
        to_owned(size());
    }

public:

    const_iterator begin() const
    {
        return m_base.m_data;
    }

    const_iterator end() const
    {
        return m_base.m_data + m_base.m_size;
    }

    const_iterator cbegin() const
    {
        return begin();
    }

    const_iterator cend() const
    {
        return end();
    }

    const_reference operator[](size_type index) const
    {
        assert(index < size());
        return m_base.m_data[index];
    }

public:

    iterator begin()
    {
        if (m_base.is_borrowed()) {
            to_owned();
        }
        return m_base.m_data;
    }

    iterator end()
    {
        if (m_base.is_borrowed()) {
            to_owned();
        }
        return m_base.m_data + m_base.m_size;
    }

    reference operator[](size_type index)
    {
        assert (index < size());
        if (m_base.is_borrowed()) {
            to_owned();
        }
        return m_base.m_data[index];
    }

public:

    constexpr size_type size() const
    {
        return m_base.m_size;
    }

    constexpr bool empty() const
    {
        return size() == 0;
    }

    constexpr vec_type type() const
    {
        return m_base.m_type;
    }

private:

    void resize_up(size_type sz)
    {
        assert(sz > size());
        if (!m_base.is_owned()) {
            to_owned(sz);
        } else {
            resize_owned(sz);
        }
    }

    void resize_up(size_type sz, const_reference val)
    {
        assert(sz > size());

        if (!m_base.is_owned()) {
            to_owned(sz, val);
        } else {
            resize_owned(sz, val);
        }

    }

    void resize_down(size_type sz)
    {
        assert(sz < size());
        if (!m_base.is_owned()) {
            to_owned(sz);
        } else {
            resize_owned(sz);
        }
    }

public:

    void resize(size_type sz)
    {
        size_type csz = size();
        if (sz > csz) {
            resize_up(sz, value_type());
        } else if (sz < csz) {
            resize_down(sz);
        }
        assert(size() == sz);
    }

    void resize(size_type sz, const_reference val)
    {
        size_type csz = size();
        if (sz > csz) {
            resize_up(sz, val);
        } else if (sz < csz) {
            resize_down(sz);
        }
        assert(size() == sz);
    }


public:

    void reserve(size_type sz)
    {
        assert (sz > size());

        if (!m_base.is_owned()) {
            to_owned(sz);
        } else {
            resize_owned(sz);
        }
        //reserve_fill(old_size, is_pod_t());
    }

    void clear()
    {
        resize(0);
    }

    void swap(dense_storage &other)
    {
        std::swap(m_base, other.m_base);
    }


private:

    void copy_extend_impl(pointer start, const_pointer begin, const_pointer end, std::true_type)
    {
        std::uninitialized_copy(begin, end, start);
    }

    void copy_extend_impl(pointer start, const_pointer begin, const_pointer end, std::false_type)
    {
        std::copy(begin, end, start);
    }

    void move_extend_impl(pointer start, pointer begin, pointer end, std::true_type)
    {
        std::uninitialized_copy(std::make_move_iterator(begin), std::make_move_iterator(end), start);
    }

    void move_extend_impl(pointer start, pointer begin, pointer end, std::false_type)
    {
        std::copy(std::make_move_iterator(begin), std::make_move_iterator(end), start);
    }


public:

    void copy_extend(const_pointer start_ptr, const_pointer end_ptr)
    {
        if (start_ptr == end_ptr) {
            return;
        }

        size_type old_size = size();
        size_type new_size = old_size + static_cast<size_type>(end_ptr - start_ptr);

        base_type new_base{new_size};

        if (m_base.is_owned()) {
            std::uninitialized_copy(
                    std::make_move_iterator(begin()),
                    std::make_move_iterator(end()),
                    new_base.m_data
                    );
        } else {
            std::uninitialized_copy(
                    begin(),
                    end(),
                    new_base.m_data
                    );
        }
        std::uninitialized_copy(start_ptr, end_ptr, new_base.m_data + m_base.m_size);
        m_base = std::move(new_base);
    }

    void move_extend(pointer start_ptr, pointer end_ptr)
    {
        if (start_ptr == end_ptr) {
            return;
        }

        size_type old_size = size();
        size_type new_size = old_size + static_cast<size_type>(end_ptr - start_ptr);

        base_type new_base{new_size};

        if (m_base.is_owned()) {
            std::uninitialized_copy(
                    std::make_move_iterator(begin()),
                    std::make_move_iterator(end()),
                    new_base.m_data
                    );
        } else {
            std::uninitialized_copy(
                    begin(),
                    end(),
                    new_base.m_data
                    );
        }
        std::uninitialized_copy(
                std::make_move_iterator(start_ptr),
                std::make_move_iterator(end_ptr),
                new_base.m_data+m_base.m_size);
        m_base = std::move(new_base);
    }

public:

    bool operator==(dense_storage const &other) const
    {
        if (size() != other.size()) {
            return false;
        }

        for (size_type i = 0; i < size(); ++i) {
            if (m_base.m_data[i] != other.m_base.m_data[i]) {
                return false;
            }
        }

        return true;
    }

    friend std::ostream &operator<<(std::ostream &os, dense_storage const &arg)
    {
        os << '{';
        for (size_type i = 0; i < arg.size(); ++i) {
            os << ' ' << arg[i];
        }
        os << " }";
        return os;
    }


};


} // namespace vectors
} // namespace alg
#endif //LIBALGEBRA_DENSE_STORAGE_H
