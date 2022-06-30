//
// Created by sam on 17/02/2021.
//

#ifndef LIBALGEBRA_UTILS_ITERATORS_H
#define LIBALGEBRA_UTILS_ITERATORS_H

#include <iterator>
#include <utility>
#include <type_traits>

#include "libalgebra/utils/meta.h"

namespace alg {
namespace utils {
namespace iterators {

template<typename Vector, typename Iterator>
inline typename std::enable_if<std::is_same<typename std::remove_cv<
            typename std::iterator_traits<Iterator>::value_type::first_type>::type,
            typename Vector::KEY>::value &&
       std::is_same<typename std::remove_cv<
            typename std::iterator_traits<Iterator>::value_type::second_type>::type,
            typename Vector::SCALAR>::value,
        typename std::iterator_traits<Iterator>::value_type::first_type
>::type
key(Iterator& it)
{
    return it->first;
}

template<typename Vector, typename Iterator>
inline typename std::enable_if<std::is_same<typename std::remove_cv<
        typename std::iterator_traits<Iterator>::value_type::first_type>::type,
                                            typename Vector::KEY>::value
                                       && std::is_same<typename std::remove_cv<
                                               typename std::iterator_traits<Iterator>::value_type::second_type>::type,
                                                       typename Vector::SCALAR>::value,
                               typename copy_constness<typename std::iterator_traits<Iterator>::value_type, typename Vector::SCALAR>::type&
           >::type
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
inline typename ValueType::value_type value(alg::vectors::iterators::vector_iterator<ValueType>& it)
{
    return it->value();
}

}// namespace iterators
}// namespace utils
}// namespace alg

#endif// LIBALGEBRA_UTILS_ITERATORS_H
