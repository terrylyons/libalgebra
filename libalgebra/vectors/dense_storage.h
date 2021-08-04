//
// Created by sam on 04/08/2021.
//

#ifndef LIBALGEBRA_DENSE_STORAGE_H
#define LIBALGEBRA_DENSE_STORAGE_H

#include <memory>

namespace alg {
namespace vectors {


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
template <typename S> class dense_storage
{
public:

    enum vec_type
    {
        owned, borrowed_mut, borrowed
    };

    using size_type = DIMN;
    using value_type = S;
    using reference = S &;
    using const_reference = S const &;
    using pointer = S *;
    using const_pointer = S const *;

    using iterator = S *;
    using const_iterator = S const *;


private:
    vec_type m_type;
    pointer m_data;
    size_type m_size;
    /* size_type m_capacity; */

public:


    dense_storage() : m_data{nullptr}, m_type{owned}, m_size{0} {}

    explicit dense_storage(size_type size) : m_data{nullptr}, m_type{owned}, m_size{size}
    {
        if (size > 0) {
            m_data = alloc_new(size);
        }
    }

    dense_storage(dense_storage const &other) : m_type{owned}, m_data{nullptr}, m_size{}
    {
        if (other.m_data) {
            alloc_new_and_copy(other.m_data, other.m_size, other.m_size);
        }
    }

    dense_storage(pointer ptr, size_type size) : m_data{ptr}, m_size{size}, m_type{borrowed_mut} {}

    dense_storage(const_pointer ptr, size_type size) : m_data{const_cast<pointer>(ptr)}, m_size{size},
                                                       m_type{borrowed} {}

    dense_storage(const_pointer begin, const_pointer end)
        : m_data{const_cast<pointer>(begin)}, m_size{static_cast<size_type>(end-begin)}, m_type{borrowed_mut}
    {}

    dense_storage(size_type offset, const_pointer start, const_pointer end)
        : m_data{}, m_size{0}, m_type{owned}
    {
        size_type size = offset + static_cast<size_type>(end - start);
        alloc_new(size);
        fill_range_default_construct(m_data, m_data+offset);
        std::uninitialized_copy(start, end, m_data+offset);
    }

    dense_storage(size_type offset, pointer start, pointer end)
        : m_data{}, m_size{0}, m_type{owned}
    {
        size_type size = offset + static_cast<size_type>(end - start);
        alloc_new(size);
        fill_range_default_construct(m_data, m_data+offset);
        std::uninitialized_copy(
                std::make_move_iterator(start),
                std::make_move_iterator(end),
                m_data+offset
                );
    }

    ~dense_storage()
    {
        if (m_type == owned && m_data) {
            dealloc(m_data, m_size);
        }
        m_data = nullptr;
        m_size = 0;
        m_type = owned;
    }

    dense_storage &operator=(dense_storage const &other)
    {
        m_type = owned;
        if (other.m_data) {
            alloc_new_and_copy(other.m_data, other.m_size, other.m_size);
        }
        return *this;
    }

    dense_storage& operator=(dense_storage&& other) noexcept
    {
        if (m_type == owned) {
            destroy(m_data, m_size);
        }
        m_data = other.m_data;
        m_size = other.m_size;
        m_type = other.m_type;
    }


private:

    void alloc_new(size_type size)
    {
        if (size > 0) {
            pointer old_data = m_data;
            try {
                m_data = new S[size];
            } catch (...) {
                m_data = old_data;
                throw;
            }
            m_size = size;
        } else {
            m_data = nullptr;
            m_size = 0;
        }
    }

    void realloc_with_copy(const_pointer old_data, size_type old_size, size_type new_size)
    {
        alloc_new(new_size);
        if (old_data) {
            std::uninitialized_copy(old_data, old_data + std::min(new_size, old_size), m_data);
        }
    }

    void realloc_with_move(pointer old_data, size_type old_size, size_type new_size)
    {
        alloc_new(new_size);
        if (old_data) {
            std::uninitialized_copy(std::make_move_iterator(old_data),
                                    std::make_move_iterator(old_data + std::min(old_size, new_size)), m_data);
        }
    }

    static void fill_range_default_construct(pointer range_start, pointer range_end)
    {
        for (pointer p = range_start; p != range_end; ++p) {
            new(p) S();
        }
    }

    static void fill_range_copy(pointer range_start, pointer range_end, const_reference val)
    {
        std::uninitialized_fill(range_start, range_end, val);
    }

    void alloc_new_and_copy(const_pointer old_data, size_type old_size, size_type new_size)
    {
        realloc_with_copy(old_data, old_size, new_size);
        if (new_size > old_size) {
            fill_range_default_construct(m_data + old_size, m_data + new_size);
        }
    }

    void alloc_new_and_copy(const_pointer old_data, size_type old_size, size_type new_size, const_reference val)
    {
        realloc_with_copy(old_data, old_size, new_size);
        if (new_size > old_size) {
            fill_range_copy(m_data + old_size, m_data + new_size, val);
        }
    }

    void alloc_new_and_move(pointer old_data, size_type old_size, size_type new_size)
    {
        realloc_with_move(old_data, old_size, new_size);

        if (new_size > old_size) {
            fill_range_default_construct(m_data + old_size, m_data + new_size);
        }

    }

    void alloc_new_and_move(pointer old_data, size_type old_size, size_type new_size, const_reference val)
    {
        realloc_with_move(old_data, old_size, new_size);
        if (new_size > old_size) {
            fill_range_copy(m_data + old_size, m_data + new_size, val);
        }
    }


    static void destroy_range(pointer, pointer, std::true_type) {}

    static void destroy_range(pointer start, pointer end, std::false_type)
    {
        while (end != start) {
            --end;
            end->~S();
        }
    }

    void destroy(pointer data, size_type size)
    {
        //destroy_range(data, data+size, std::is_trivially_destructible<S>());
        delete[] data;
    }

    void dealloc(pointer data, size_type size)
    {
        assert(m_type == owned);
        assert(m_data != nullptr);
        destroy(m_data, m_size);
        m_data = nullptr;
        m_size = 0;
        m_type = owned;
        //alloc_traits::deallocate(m_alloc, data, size);
    }

    void to_owned(size_type alloc_size)
    {
        assert(m_type != owned);
        m_type = owned;
        alloc_new_and_copy(m_data, m_size, alloc_size);
    }

    void to_owned(size_type alloc_size, const_reference val)
    {
        assert(m_type != owned);
        m_type = owned;
        alloc_new_and_copy(m_data, m_size, alloc_size, val);
    }

    void resize_owned(size_type size)
    {
        assert(m_type == owned);
        assert(m_size != size);

        pointer old_data = m_data;
        DIMN const old_size = m_size;
        alloc_new_and_move(old_data, old_size, size);
        assert(m_data != old_data);
        destroy(old_data, old_size);
    }

    void resize_owned(size_type size, const_reference val)
    {
        assert(m_type == owned);
        assert(m_size != size);

        pointer old_data = m_data;
        DIMN const old_size = m_size;
        alloc_new_and_move(old_data, old_size, size, val);
        assert(m_data != old_data);
        destroy(old_data, old_size);
    }

    void to_owned()
    {
        to_owned(size());
    }

public:

    const_iterator begin() const
    {
        return m_data;
    }

    const_iterator end() const
    {
        return m_data + m_size;
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
        return m_data[index];
    }

public:

    iterator begin()
    {
        if (m_type == borrowed) {
            to_owned();
        }
        return m_data;
    }

    iterator end()
    {
        if (m_type == borrowed) {
            to_owned();
        }
        return m_data + m_size;
    }

    reference operator[](size_type index)
    {
        assert (index < size());
        if (m_type == borrowed) {
            to_owned();
        }
        return m_data[index];
    }

public:

    constexpr size_type size() const
    {
        return m_size;
    }

    constexpr bool empty() const
    {
        return size() == 0;
    }

    constexpr vec_type type() const
    {
        return m_type;
    }

private:

    void resize_up(size_type size)
    {
        assert(size > m_size);
        if (m_type != owned) {
            to_owned(size);
        } else {
            resize_owned(size);
        }
    }

    void resize_up(size_type size, const_reference val)
    {
        assert(size > m_size);
        if (m_type != owned) {
            to_owned(size, val);
        } else {
            resize_owned(size, val);
        }
    }

    void resize_down(size_type size)
    {
        assert(size < m_size);
        if (m_type != owned) {
            alloc_new_and_copy(m_data, m_size, size);
        } else {
            resize_owned(size);
        }
    }

public:

    void resize(size_type size)
    {
        if (size > m_size) {
            resize_up(size);
        } else if (size < m_size) {
            resize_down(size);
        }
        assert(m_size == size);
    }

    void resize(size_type size, const_reference val)
    {
        if (size > m_size) {
            resize_up(size, val);
        } else if (size < m_size) {
            resize_down(size);
        }
        assert(m_size == size);
    }

    void clear()
    {
        if (m_type == owned) {
            dealloc(m_data, m_size);
        }
        m_data = nullptr;
        m_size = 0;
        m_type = owned;
    }

    void swap(dense_storage &other)
    {
        std::swap(m_data, other.m_data);
        std::swap(m_size, other.m_size);
        std::swap(m_type, other.m_type);
    }


public:

    void copy_extend(const_pointer start, const_pointer end)
    {
        if (start != end) {
            size_type old_size = m_size;
            size_type new_size = m_size + static_cast<size_type>(end - start);
            realloc_with_move(m_data, m_size, new_size);
            assert(m_size == new_size);
            std::uninitialized_copy(start, end, m_data + old_size);
        }
    }

    void move_extend(pointer start, pointer end)
    {
        if (start != end) {
            size_type old_size = m_size;
            size_type new_size = m_size + static_cast<size_type>(end - start);
            realloc_with_move(m_data, m_size, new_size);
            assert(m_size == new_size);
            std::uninitialized_copy(
                    std::make_move_iterator(start),
                    std::make_move_iterator(end),
                    m_data + old_size
                    );
        }
    }


};


}; // namespace vectors
}; // namespace alg
#endif //LIBALGEBRA_DENSE_STORAGE_H
