//
// Created by sam on 18/02/2021.
//

#ifndef LIBALGEBRAUNITTESTS_COMPAT_H
#define LIBALGEBRAUNITTESTS_COMPAT_H

#include <iterator>
#include <utility>

#include <boost/type_traits.hpp>
#include <boost/mpl/logical.hpp>

#include <libalgebra/detail/meta.h>
#include <libalgebra/iterators.h>



namespace iter {

#define ENABLE_IF_PAIR(OTYPE, ITER, KT, VT)                                     \
    typename boost::enable_if<                                                  \
        typename boost::mpl::and_<                                              \
            typename boost::is_same<                                            \
                typename boost::remove_cv<                                      \
                    typename std::iterator_traits<Iterator>::value_type::first_type                                 \
                    >::type,                                                    \
                KT                                                              \
            >::type,                                                            \
            typename boost::is_same<                                            \
                typename boost::remove_cv<                                      \
                    typename std::iterator_traits<Iterator>::value_type::second_type                                \
                >::type,                                                        \
                VT                                                              \
            >::type                                                             \
        >::type,                                                                \
        OTYPE                                                                   \
    >::type


template <typename Vector, typename Iterator>
ENABLE_IF_PAIR(typename Vector::KEY, Iterator, typename Vector::KEY, typename Vector::SCALAR)
key(Iterator& it)
{ return it->first; }

template <typename Vector, typename Iterator>
typename alg::utils::copy_constness<
        typename std::iterator_traits<Iterator>::reference,
        ENABLE_IF_PAIR(typename Vector::SCALAR, Iterator, typename Vector::KEY, typename Vector::SCALAR)
>::type
value(Iterator& it)
{ return it->second; }



#undef ENABLE_IF_PAIR

template <typename Vector, typename ValueType>
typename ValueType::key_type key(alg::vectors::iterators::vector_iterator<ValueType>& it)
{
    return it->key();
}

template <typename Vector, typename ValueType>
typename ValueType::value_type value(alg::vectors::iterators::vector_iterator<ValueType>& it)
{
    return it->value();
}


}


#endif //LIBALGEBRAUNITTESTS_COMPAT_H
