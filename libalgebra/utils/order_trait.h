//
// Created by sam on 30/03/2021.
//

#ifndef LIBALGEBRA_ORDER_TRAIT_H
#define LIBALGEBRA_ORDER_TRAIT_H

namespace alg {
namespace utils {

template <typename Map>
struct is_ordered {
    static const bool value = false;
};

template <typename K,
        typename V>
struct is_ordered<std::map<K,
                           V> > {
    static const bool value = true;
};

}
}


#endif //LIBALGEBRA_ORDER_TRAIT_H
