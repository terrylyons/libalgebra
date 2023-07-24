//
// Created by sam on 04/08/2021.
//

#ifndef LIBALGEBRA_DENSE_STORAGE_H
#define LIBALGEBRA_DENSE_STORAGE_H

#include <cassert>
#include <initializer_list>
#include <iostream>
#include <memory>

#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
#include <boost/serialization/array.hpp>
#endif
#include "libalgebra/implementation_types.h"

namespace alg {
namespace vectors {

namespace dtl {

/**
 * @brief Base storage for dense vectors.
 *
 * This is the base storage type for dense vectors. It handles allocation and
 * deallocating space for a dense_storage object.
 *
 * @tparam S Scalar type to be held in the vector.
 * @tparam Alloc Allocator type to use for allocation.
 */
template<typename S, typename Alloc>
struct dense_storage_base {
    using allocator_type = Alloc;
    using alloc_traits = std::allocator_traits<allocator_type>;
    using value_type = S;
    using pointer = S*;
    using const_pointer = S const*;
    using size_type = typename allocator_type::size_type;

    /**
     * The vector can be owned where the container is responsible for the
     * data it points to; borrowed where the container points to some data
     * in a const way; or borrowed_mut where data is borrowed in a non const
     * way.
     */
    enum vec_type
    {
        owned,
        borrowed_mut,
        borrowed
    };

    allocator_type m_alloc;
    pointer m_data;
    size_type m_size;
    vec_type m_type;

    /// Create new storage (default initialised) with size
    explicit dense_storage_base(size_type sz = 0)
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

    /// Create a new mutably borrowed vector from data
    dense_storage_base(pointer begin, pointer end)
        : m_alloc{}, m_data{begin}, m_size{static_cast<size_type>(end - begin)}, m_type{borrowed_mut}
    {}

    /// Create a new borrowed data from data
    dense_storage_base(const_pointer begin, const_pointer end)
        : m_alloc{},
          m_data{const_cast<pointer>(begin)},
          m_size{static_cast<size_type>(end - begin)},
          m_type{borrowed}
    {}

    // No copy operations
    dense_storage_base(dense_storage_base const&) = delete;
    dense_storage_base& operator=(dense_storage_base const&) = delete;

    /// Move constructor
    dense_storage_base(dense_storage_base&& other) noexcept
        : m_alloc{other.m_alloc}, m_data{other.m_data}, m_size{other.m_size}, m_type{other.m_type}
    {
        other.m_data = nullptr;
        other.m_size = 0;
        other.m_type = owned;
    }

    /// Move assignment
    dense_storage_base& operator=(dense_storage_base&& other) noexcept
    {
        if (this != &other) {
            if (m_type == owned) {
                alloc_traits::deallocate(m_alloc, m_data, m_size);
            }
            m_alloc = std::move(other.m_alloc);
            m_data = other.m_data;
            m_type = other.m_type;
            m_size = other.m_size;

            other.m_data = nullptr;
            other.m_size = 0;
        }

        return *this;
    }

    /// Test if storage is owned
    bool is_owned() const
    {
        return m_type == owned;
    }

    /// Test if storage is borrowed
    bool is_borrowed() const
    {
        return m_type == borrowed;
    }

    /// Test if storage is borrowed mutably
    bool is_borrowed_mut() const
    {
        return m_type == borrowed_mut;
    }

    template<typename... Args>
    value_type& emplace(size_type idx, Args&&... args) noexcept(noexcept(alloc_traits::construct(m_alloc, m_data + idx, std::forward<Args>(args)...)))
    {
        assert(idx < m_size);
        alloc_traits::construct(m_alloc, m_data + idx, std::forward<Args>(args)...);
        return m_data[idx];
    }
};

}// namespace dtl

/**
 * @brief Storage class for dense vectors.
 *
 * This class will handle the storage for the dense vector type. For the time being, this will simply
 * hold a pointer to the data, the size of the buffer and the type of the vector which can be either
 * `owned`, `borrowed_mut`, `borrowed`. This storage type will implement a "copy on resize" for `borrowed`
 * and `borrowed_mut` types, and "copy on modify" for "borrowed" types.
 *
 * Fundamentally, this is a clone of the C++ standard libary vector type, with an internal state which
 * describes whether the data is owned or borrowed. However, there are some key differences. Most importantly
 * the reserve member function allocates new memory and resizes the vector but it does not instantiate the
 * elements within the new memory (unless they are copied/moved from old data). This is useful for situations
 * where the elements are to be immediately overwritten.
 *
 * @tparam S Scalar type to store.
 * @tparam Alloc Allocator to use for allocating and deallocating vectors. Default is std::allocator<S>.
 */
template<typename S, typename Alloc = std::allocator<S>>
class dense_storage
{
    using base_type = dtl::dense_storage_base<S, Alloc>;

public:
    using allocator_type = typename base_type::allocator_type;
    using alloc_traits = typename base_type::alloc_traits;

    using size_type = typename base_type::size_type;


    using value_type = S;
    using reference = S&;
    using const_reference = S const&;
    using pointer = S*;
    using const_pointer = S const*;

    using iterator = S*;
    using const_iterator = S const*;

    using vec_type = typename base_type::vec_type;

private:
    base_type m_base;

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
        std::is_trivially_destructible<S> tag;
        destroy_range(start, end, tag);
    }

    void fill_range_default_construct(pointer start, pointer end)
    {
        std::uninitialized_fill(start, end, value_type());
    }

public:
    /**
     * @brief Constructor for blank data with given size
     *
     * Creates a new storage with the given size. Elements are default initialised.
     *
     * @param size Size of buffer to allocate
     */
    explicit dense_storage(size_type size = 0)
        : m_base(size)
    {
        fill_range_default_construct(m_base.m_data, m_base.m_data + m_base.m_size);
    }

    /**
     * @brief Copy constructor
     *
     * Makes a new copy of other. The resulting storage owns its data even if other does not.
     *
     * @param other Storage to copy data from
     */
    dense_storage(dense_storage const& other)
        : m_base(other.size())
    {
        if (m_base.m_data != nullptr) {
            assert(other.begin() != nullptr && other.end() != nullptr);
            std::uninitialized_copy(other.begin(), other.end(), m_base.m_data);
        }
    }

    /**
     * @brief Constructor from mutable pointer to existing data
     *
     * Creates a new mutably borrowed storage pointing to the range provided.
     *
     * @param ptr Start of data range to borrow
     * @param size Size of data range to borrow
     */
    dense_storage(pointer ptr, size_type size)
        : m_base{ptr, ptr + size}
    {
    }

    /**
     * @brief Constructor from const pointer to existing data
     *
     * Creates a new borrowed storage pointing to the range provided.
     *
     * @param ptr Start of data range to borrow
     * @param size Size of data range to borrow
     */
    dense_storage(const_pointer ptr, size_type size)
        : m_base{ptr, ptr + size}
    {
    }

    /**
     * @brief Constructor from mutable pointer to existing data range
     *
     * Construct a new mutably borrowed storage from pointers to range [begin, end).
     *
     * @param begin Pointer to beginning of range to borrow
     * @param end Pointer to one past end of the range to borrow.
     */
    dense_storage(pointer begin, pointer end)
        : m_base{begin, end}
    {
    }

    /**
     * @brief Constructor from const pointer to existing data range
     *
     * Construct a new borrowed storage from pointers to range [begin, end).
     *
     * @param begin Pointer to beginning of range to borrow
     * @param end Pointer to one past end of the range to borrow.
     */
    dense_storage(const_pointer begin, const_pointer end)
        : m_base{begin, end}
    {
    }

    /**
     * @brief Constructor for new owned storage with offset followed by existing data
     *
     * Create a new owned storage of size offset + (end - start) where the first offset
     * elements are default initialised and the remaining buffer is filled with the
     * values from [start, end).
     *
     * @param offset Size of offset to prepend to storage
     * @param start start of data range to copy data from
     * @param end One past end of data range to copy data range from
     */
    dense_storage(size_type offset, const_pointer start, const_pointer end)
        : m_base{offset + static_cast<size_type>(end - start)}
    {
        fill_range_default_construct(m_base.m_data, m_base.m_data + offset);
        std::uninitialized_copy(start, end, m_base.m_data + offset);
    }

    /**
     * @brief Constructor for new owned storage with offset followed by existing data
     *
     * Create a new owned storage of size offset + (end - start) where the first offset
     * elements are default initialised and the remaining buffer is filled with the
     * values from [start, end).
     *
     * @param offset Size of offset to prepend to storage
     * @param start start of data range to copy data from
     * @param end One past end of data range to copy data range from
     */
    dense_storage(size_type offset, pointer start, pointer end)
        : m_base{offset + static_cast<size_type>(end - start)}
    {
        fill_range_default_construct(m_base.m_data, m_base.m_data + offset);
        std::uninitialized_copy(std::make_move_iterator(start), std::make_move_iterator(end), m_base.m_data + offset);
    }

    /**
     * @brief Initializer list constructor
     *
     */
    dense_storage(std::initializer_list<S> args) : m_base(args.size())
    {
        size_type i = 0;
        for (auto v : args) {
            m_base.emplace(i++, v);
        }
    }

    ~dense_storage()
    {
        if (m_base.is_owned()) {
            destroy_range(m_base.m_data, m_base.m_data + m_base.m_size);
        }
    }

    /// Copy constructor - the new storage owns its data even if other borrows data.
    dense_storage& operator=(dense_storage const& other)
    {
        //        dense_storage tmp(other);
        //        this->swap(tmp);
        if (m_base.is_owned()) {
            destroy_range(m_base.m_data, m_base.m_data + m_base.m_size);
        }
        m_base = base_type(other.m_base.m_size);
        std::uninitialized_copy(other.m_base.m_data, other.m_base.m_data + other.m_base.m_size, m_base.m_data);
        return *this;
    }

    /// Move constructor
    dense_storage& operator=(dense_storage&& other) noexcept
    {
        if (m_base.is_owned()) {
            destroy_range(m_base.m_data, m_base.m_data + m_base.m_size);
        }
        m_base = std::move(other.m_base);
        return *this;
    }

    /**
     * @brief Copy stored data to owned data into a new buffer
     *
     * Allocates new owned storage copies a range of existing data into the new buffer.
     *
     * @param ptr Pointer to start of existing data.
     * @param sz Size of existing data buffer buffer.
     * @return New owned dense_storage containing a copy of the data
     */
    static dense_storage make_owned(const_pointer ptr, size_type sz)
    {
        dense_storage result(ptr, sz);
        result.to_owned();
        return result;
    }

    operator bool() const noexcept
    {
        return m_base.m_size != 0;
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

        if (mid > 0) {
            std::uninitialized_copy(m_base.m_data, m_base.m_data + mid, new_base.m_data);
        }
        if (alloc_size > mid) {
            maybe_fill(new_base.m_data + mid, new_base.m_data + new_base.m_size, value_type());
        }

        m_base = std::move(new_base);
    }

    void to_owned(size_type alloc_size, const_reference val)
    {
        assert(!m_base.is_owned());
        base_type new_base(alloc_size);
        size_type mid = std::min(alloc_size, size());

        std::uninitialized_copy(m_base.m_data, m_base.m_data + m_base.m_size, new_base.m_data);

        if (mid < alloc_size) {
            // new larger than old
            std::uninitialized_fill(new_base.m_data + m_base.m_size, new_base.m_data + alloc_size, val);
        }

        m_base = std::move(new_base);
    }

    void resize_owned(size_type sz)
    {
        if (sz == 0) {
            m_base = base_type();
            return;
        }

        assert(m_base.is_owned());
        assert(size() != sz);

        size_type mid = std::min(sz, size());

        base_type new_base(sz);
        std::uninitialized_copy(
                std::make_move_iterator(m_base.m_data),
                std::make_move_iterator(m_base.m_data + mid),
                new_base.m_data);

        if (mid >= sz) {
            destroy_range(m_base.m_data + mid, m_base.m_data + m_base.m_size);
        }
        else {
            maybe_fill(new_base.m_data + mid, new_base.m_data + new_base.m_size, value_type());
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
                std::make_move_iterator(m_base.m_data + mid),
                new_base.m_data);

        if (mid < size()) {
            destroy_range(m_base.m_data + mid, m_base.m_data + m_base.m_size);
        }
        else {
            std::uninitialized_fill(new_base.m_data + mid, new_base.m_data + sz, val);
        }

        m_base = std::move(new_base);
    }

    void to_owned()
    {
        to_owned(size());
    }

public:
    /// Const pointer to beginning of storage
    const_iterator begin() const
    {
        return m_base.m_data;
    }

    /// Const pointer to one past end of storage
    const_iterator end() const
    {
        return m_base.m_data + m_base.m_size;
    }

    /// Const pointer to beginning of storage
    const_iterator cbegin() const
    {
        return begin();
    }

    /// Const pointer to one past end of storage
    const_iterator cend() const
    {
        return end();
    }

    /// Index access to const data
    const_reference operator[](size_type index) const
    {
        assert(index < size());
        return m_base.m_data[index];
    }

public:
    /**
     * @brief Pointer to beginning of storage
     *
     * For owned or mutably-borrowed data this simply returns a pointer to the
     * start of the storage range. For const borrowed storage, the data is first
     * copied into an owned buffer.
     *
     * @return (mutable) pointer to beginning of storage
     */
    iterator begin()
    {
        if (m_base.is_borrowed()) {
            to_owned();
        }
        return m_base.m_data;
    }

    /**
     * @brief Pointer to one past end of storage
     *
     * For owned or mutably-borrowed data this simply returns a pointer to one past
     * the end of the storage range. For const borrowed storage, the data is first
     * copied into an owned buffer.
     *
     * @return (mutable) pointer to one past the end of storage
     */
    iterator end()
    {
        if (m_base.is_borrowed()) {
            to_owned();
        }
        return m_base.m_data + m_base.m_size;
    }

    /**
     * @brief Mutable index access
     *
     * For owned or mutably-borrowed data this simply returns a pointer to the
     * element at index of the storage range. For const borrowed storage, the
     * data is first copied into an owned buffer.
     *
     * Index bounds are not checked for performance.
     *
     * @param index index in storage of element to access
     * @return mutable reference to element at index in storage.
     */
    reference operator[](size_type index)
    {
        assert(index < size());
        if (m_base.is_borrowed()) {
            to_owned();
        }
        return m_base.m_data[index];
    }

public:
    /// Get the size of storage
    constexpr size_type size() const
    {
        return m_base.m_size;
    }

    /// Check if storage is empty
    constexpr bool empty() const
    {
        return size() == 0;
    }

    /**
     * @brief Get the type of storage
     *
     * @return one of borrowed, borrowed_mut, or owned.
     */
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
        }
        else {
            resize_owned(sz);
        }
    }

    void resize_up(size_type sz, const_reference val)
    {
        assert(sz > size());

        if (!m_base.is_owned()) {
            to_owned(sz, val);
        }
        else {
            resize_owned(sz, val);
        }
    }

    void resize_down(size_type sz)
    {
        assert(sz < size());
        if (!m_base.is_owned()) {
            to_owned(sz);
        }
        else {
            resize_owned(sz);
        }
    }

public:
    /**
     * @brief Resize the storage to new size.
     *
     * Converts the storage to owned data and resizes to new size.
     * If the new size is larger, the new elements are default initialised.
     *
     * @param sz new size for the storage
     */
    void resize(size_type sz)
    {
        size_type csz = size();
        if (sz > csz) {
            resize_up(sz, value_type());
        }
        else if (sz < csz) {
            resize_down(sz);
        }
        assert(size() == sz);
    }

    /**
     * @brief Resize the storage to new size with specified fill value
     *
     * Converts the storage to owned data and resizes to new size.
     * If the new storage is larger, the new elements are initialised to val.
     *
     * @param sz new size for the storage
     * @param val value to fill new elements (if any) with
     */
    void resize(size_type sz, const_reference val)
    {
        size_type csz = size();
        if (sz > csz) {
            resize_up(sz, val);
        }
        else if (sz < csz) {
            resize_down(sz);
        }
        assert(size() == sz);
    }

public:
    /**
     * @brief Grow the storage without initialising new values
     *
     * Converts the data to owned if it wasn't already, where the new buffer
     * has the target size. The excess buffer space is not filled with any values
     * and it is assumed that you will initialise these values very shortly after
     * performing a reserve.
     *
     * @param sz new target size
     */
    void reserve(size_type sz)
    {
        if (sz <= size()) {
            return;
        }

        if (!m_base.is_owned()) {
            to_owned(sz);
        }
        else {
            resize_owned(sz);
        }
        //reserve_fill(old_size, is_pod_t());
    }

    /**
     * @brief Clear the storage
     *
     * For borrowed data, this is simply "forgetting" the data that it points to.
     * For owned data, the storage is cleared and deallocated.
     */
    void clear()
    {
        resize(0);
    }

    /// Swap this storage with another
    void swap(dense_storage& other)
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
    /**
     * @brief Extend the storage by copying data from the range [start_ptr, end_ptr)
     *
     * Make the storage owned and grow the allocated space by (end_ptr - start_ptr).
     * New space is filled by copying the data from the range [start_ptr, end_ptr).
     *
     * @param start_ptr Start of range to copy into new storage
     * @param end_ptr One past end of range to copy into new storage
     */
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
                    new_base.m_data);
        }
        else {
            std::uninitialized_copy(
                    begin(),
                    end(),
                    new_base.m_data);
        }
        std::uninitialized_copy(start_ptr, end_ptr, new_base.m_data + m_base.m_size);
        m_base = std::move(new_base);
    }

    /**
     * @brief Extend the storage by moving the data from the range [start_ptr, end_ptr)
     *
     * Make the storage owned and grow the allocated space by (end_ptr - start_ptr).
     * New space is filled by moving the data from the range [start_ptr, end_ptr).
     *
     * @param start_ptr Start of range to move into new storage
     * @param end_ptr one past end of range to move into new storage
     */
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
                    new_base.m_data);
        }
        else {
            std::uninitialized_copy(
                    begin(),
                    end(),
                    new_base.m_data);
        }
        std::uninitialized_copy(
                std::make_move_iterator(start_ptr),
                std::make_move_iterator(end_ptr),
                new_base.m_data + m_base.m_size);
        m_base = std::move(new_base);
    }

    template<typename... Args>
    reference emplace(size_type i, Args&&... args) noexcept(noexcept(value_type(args...)))
    {
        return m_base.emplace(i, std::forward<Args>(args)...);
    }

public:
    /// Equality operator
    bool operator==(dense_storage const& other) const
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

    /// Print data to stream (useful for debugging)
    friend std::ostream& operator<<(std::ostream& os, dense_storage const& arg)
    {
        os << '{';
        for (size_type i = 0; i < arg.size(); ++i) {
            os << ' ' << arg[i];
        }
        os << " }";
        return os;
    }

#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
private:
    friend class boost::serialization::access;

    BOOST_SERIALIZATION_SPLIT_MEMBER()

    template<typename Archive>
    void load(Archive& ar, const unsigned /* version */)
    {
        size_type sz = 0;
        ar >> sz;
        m_base = base_type(sz);
        if (sz > 0) {
            assert(m_base.m_data != nullptr);
            fill_range_default_construct(m_base.m_data, m_base.m_data + m_base.m_size);
            ar >> boost::serialization::make_array(m_base.m_data, m_base.m_size);
        }
    }

    template<typename Archive>
    void save(Archive& ar, const unsigned int /*version*/) const
    {
        ar << boost::serialization::make_nvp("size", m_base.m_size);
        ar << boost::serialization::make_array(m_base.m_data, m_base.m_size);
    }
#endif
};

}// namespace vectors
}// namespace alg
#endif//LIBALGEBRA_DENSE_STORAGE_H
