//
// Created by sam on 17/02/2021.
//

#ifndef LIBALGEBRA_VECTORS_ITERATORS_H
#define LIBALGEBRA_VECTORS_ITERATORS_H

#include <iterator>

namespace alg {
namespace vectors {
namespace iterators {

/*
 * The value type of a vector iterator should have the following format.
 *
 * (Template parameters need not be actual template parameters, but might be
 * implicitly defined in the vector)
 * template <typename Vector, typename Iterator, typename Key, typename Value>
 * class vector_iterator_item
 * {
 * public:
 *      // Declare friend of vector iterator so it has access to
 *      // the internals
 *      friend class vector_iterator<vector_iterator_item>;
 *
 *      vector_iterator_item(); // Default constructor
 *      vector_iterator_item(const vector_iterator_item&); // Copy constructor
 *      vector_iterator_item((const) Vector&, Iterator it); // Vector/iterator
 * initialisation
 *
 *      // type definitions of output types
 *      typedef Key key_type;
 *      typedef Value value_type;
 *
 *      key_type key();
 *      value_type& value();
 *
 * private:
 *      Iterator m_iterator;
 *
 * private:
 *
 *      // Compare the iterators, return True on equal
 *      bool compare_iterator(const vector_iterator_item& other);
 *
 *      // Const iterator
 *      void advance_iterator();
 *
 * };
 *
 *
 */

/**
 * @brief Iterator implementation for vector types
 *
 * This is a generic wrapper that provides a unified iterator object for vector types.
 * The vector must define a ValueType, which has member functions to get the key and value
 * of the corresponding element in the vector. This class then wraps the value type
 * and implements all the standard iterator methods.
 *
 * @tparam ValueType Item type of the iterator.
 */
template<typename ValueType>
class vector_iterator
{
public:
    typedef ValueType value_type;
    typedef const value_type& reference;
    typedef const value_type* pointer;
    typedef std::ptrdiff_t difference_type;
    typedef std::forward_iterator_tag iterator_category;

private:
    value_type m_value;

public:
    // Constructors

    /// Default constructor
    vector_iterator()
        : m_value()
    {}

    /// Vector/iterator constructor
    template<typename Vector, typename Iterator>
    vector_iterator(Vector& vect, Iterator it)
        : m_value(vect, it)
    {}

    /*
    vector_iterator& operator=(const vector_iterator& other)
    {
        m_value = other.m_value;
        return *this;
    }

     */
public:
    // Iterator advance methods

    /// Prefix increment
    vector_iterator& operator++()
    {
        m_value.advance();
        return *this;
    }

    /// Postfix increment
    vector_iterator operator++(int)
    {
        vector_iterator new_it(*this);
        operator++();
        return new_it;
    }

public:
    // Iterator access

    reference operator*() const
    {
        return m_value;
    }

    pointer operator->() const
    {
        return &m_value;
    }

public:
    // Comparison operators

    bool operator==(const vector_iterator& other) const
    {
        return m_value.compare_iterators(other.m_value);
    }

    bool operator!=(const vector_iterator& other) const
    {
        return !operator==(other);
    }
};

}// namespace iterators
}// namespace vectors
}// namespace alg

#endif// LIBALGEBRA_VECTORS_ITERATORS_H
