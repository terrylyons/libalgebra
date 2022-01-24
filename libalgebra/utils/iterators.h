//
// Created by sam on 17/02/2021.
//

#ifndef LIBALGEBRA_UTILS_ITERATORS_H
#define LIBALGEBRA_UTILS_ITERATORS_H

#include <iterator>
#include <utility>

#include <boost/mpl/logical.hpp>
#include <boost/type_traits.hpp>

#include "libalgebra/utils/meta.h"

namespace alg {
namespace utils {
namespace iterators {

template<typename Vector, typename Iterator>
inline typename boost::enable_if<typename boost::mpl::and_<
                                         boost::is_same<typename boost::remove_cv<typename std::iterator_traits<Iterator>::value_type::first_type>::type,
                                                        typename Vector::KEY>,
                                         boost::is_same<
                                                 typename boost::remove_cv<typename std::iterator_traits<Iterator>::value_type::second_type>::type,
                                                 typename Vector::SCALAR>>::type,
                                 typename std::iterator_traits<Iterator>::value_type::first_type>::type
key(Iterator& it)
{
    return it->first;
}

template<typename Vector, typename Iterator, typename ValueType = typename std::iterator_traits<Iterator>::value_type,
         typename Reference = typename std::iterator_traits<Iterator>::reference>
inline typename boost::enable_if<
        typename boost::mpl::and_<boost::is_same<typename ValueType::first_type, const typename Vector::KEY>,
                                  boost::is_same<typename ValueType::second_type, typename Vector::SCALAR>>::type,
        typename copy_constness<Reference, typename ValueType::second_type>::type&>::type
value(Iterator& it)
{
    return it->second;
}

template<typename Vector, typename ValueType>
inline typename ValueType::key_type key(alg::vectors::iterators::vector_iterator<ValueType>& it)
{
    return it->key();
}

template<typename Vector, typename ValueType>
inline typename ValueType::value_type& value(alg::vectors::iterators::vector_iterator<ValueType>& it)
{
    return it->value();
}

}// namespace iterators
}// namespace utils
}// namespace alg

#endif// LIBALGEBRA_UTILS_ITERATORS_H
