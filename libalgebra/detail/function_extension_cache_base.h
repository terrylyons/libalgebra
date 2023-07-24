//
// Created by sam on 20/01/2022.
//

#ifndef LIBALGEBRA_FUNCTION_EXTENSION_CACHE_BASE_H
#define LIBALGEBRA_FUNCTION_EXTENSION_CACHE_BASE_H

#include <map>
#include <mutex>

namespace alg {

template<
        typename KeyType,
        typename OutputType,
        template<typename, typename, typename...> class MapType = std::map,
        typename... Args>
class function_extension_cache_base
{
    using map_type = MapType<KeyType, OutputType, Args...>;

    mutable map_type cache;

public:
    template<typename Fn>
    void apply_with_caching(const KeyType& key, Fn&& fn) const
    {
        static std::recursive_mutex table_lock;

        std::lock_guard<std::recursive_mutex> access(table_lock);

        typename map_type::iterator it = cache.find(key);
        if (it != cache.end()) {
            return it->second;
        }

        return cache[key] = fn(key);
    }
};

}// namespace alg

#endif//LIBALGEBRA_FUNCTION_EXTENSION_CACHE_BASE_H
