//
// Created by sam on 22/11/2021.
//

#ifndef LIBALGEBRA_CACHING_TAGS_H
#define LIBALGEBRA_CACHING_TAGS_H

namespace alg {

// Tag types for extension functions caching

/**
 * @brief No caching tag
 *
 * Using this tag indicates that no caching should be used
 * and all values should be computed recursively to leaves
 * in all computations. This is most suitable for functions
 * that are computed rarely and are not expensive.
 */
struct no_caching_tag {
};

/**
 * @brief Lazy caching tag
 *
 * Using this tag indicates that values should be cached on
 * first computation conditionally based on a predicate.
 * The predicate should be a function object that evaluates
 * true on keys that should be added to the cache and false
 * if the value should be computed dynamically every time.
 *
 * Passing a void predicate will enable lazy caching for
 * all values.
 *
 * @tparam Predicate Predicate function object
 */
template<typename Predicate>
struct lazy_cache_tag {
    Predicate predicate;
};

template<>
struct lazy_cache_tag<void> {
};

/**
 * @brief Lookup table caching tag
 *
 * Using this tag indicates that values should be computed
 * and stored in a lookup table on first call. The tag takes a
 * predicate function object, like the lazy caching tag, which
 * determines the keys that should have their values written in
 * the table. Note that the table construction will terminate at
 * the first key which fails the predicate.
 *
 * @tparam Predicate Predicate function object
 */
template<typename Predicate>
struct lookup_table_tag {
    Predicate predicate;
};


struct lazy_cache_on_object_tag {
};

}// namespace alg

#endif//LIBALGEBRA_CACHING_TAGS_H
