//
// Created by sam on 17/02/2021.
//

#ifndef LIBALGEBRA_ITERATORS_H
#define LIBALGEBRA_ITERATORS_H

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
 *      template <typename Basis, typename Coeff>
 *      friend class vector_iterator<Basis, Coeff, vector_iterator_item>;
 *
 *      vector_iterator_item(); // Default constructor
 *      vector_iterator_item(const vector_iterator_item&); // Copy constructor
 *      vector_iterator_item((const) Vector&, Iterator it); // Vector/iterator initialisation
 *
 *      Key key();
 *      Value& value();
 *
 * private:
 *      Iterator m_iterator;
 *
 * private:
 *
 *      // Compare the iterators, return True on equal
 *      bool compare_iterator(const vector_iterator_item& other);
 *
 * };
 *
 *
 */

template <typename Basis, typename Coeff, typename ValueType>
class vector_iterator
{
public:

    typedef ValueType value_type;
    typedef value_type& reference;
    typedef value_type* pointer;
    typedef std::ptrdiff_t difference_type;
    typedef std::forward_iterator_tag iterator_category;

private:
    value_type m_value;

public:

    // Constructors

    /// Default constructor
    vector_iterator() : m_value()
    {}

    /// Copy constructor
    vector_iterator(const vector_iterator& other) : m_value(other.m_value)
    {}

    /// Vector/iterator constructor
    template <typename Vector, typename Iterator>
    vector_iterator(Vector& vect, Iterator it) : m_value(vect, it)
    {}

public:

    // Iterator advance methods

    /// Prefix increment
    vector_iterator& operator++()
    {
        m_value.m_iterator++;
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

    reference operator*()
    {
        return m_value;
    }

    pointer operator->()
    {
        return &m_value;
    }

public:

    // Comparison operators

    bool operator==(const vector_iterator& other) const
    {
        return m_value.compare_iterator(other.m_value);
    }

    bool operator!=(const vector_iterator& other) const
    {
        return !m_value.compare_iterator(other.m_value);
    }


};


} // namespace iterators
} // namespace vectors
} // namespace alg


#endif //LIBALGEBRA_ITERATORS_H
