//
// Created by user on 26/10/22.
//

#ifndef LIBALGEBRA_LIBALGEBRA_DETAIL_LEVEL_WALKERS_H_
#define LIBALGEBRA_LIBALGEBRA_DETAIL_LEVEL_WALKERS_H_

#include <utility>

namespace alg {
namespace dtl {

template <template <unsigned> class MetaFunc, unsigned MaxLevel, unsigned Level=0>
struct increasing_level_walker
{
    using next = increasing_level_walker<MetaFunc, MaxLevel, Level+1>;

    template <typename... Args>
    static void eval(Args&&... args)
    {
        MetaFunc<Level>::eval(args...);
        next::eval(std::forward<Args>(args)...);
    }
};

template <template <unsigned> class MetaFunc, unsigned MaxLevel>
struct increasing_level_walker<MetaFunc, MaxLevel, MaxLevel>
{
    template <typename... Args>
    static void eval(Args&&... args)
    {
        MetaFunc<MaxLevel>::eval(args...);
    }
};


template <template <unsigned> class MetaFunc, unsigned Level, unsigned MinLevel=0U>
struct decreasing_level_walker
{
    using next = decreasing_level_walker<MetaFunc, Level-1, MinLevel>;

    template <typename... Args>
    static void eval(Args&&... args)
    {
        MetaFunc<Level>::eval(args...);
        next::eval(std::forward<Args>(args)...);
    }
};

template <template <unsigned> class MetaFunc, unsigned MinLevel>
struct decreasing_level_walker<MetaFunc, MinLevel, MinLevel>
{
    template <typename... Args>
    static void eval(Args&&... args)
    {
        MetaFunc<MinLevel>::eval(args...);
    }
};




} // namespace dtl
} // namespace alg


#endif//LIBALGEBRA_LIBALGEBRA_DETAIL_LEVEL_WALKERS_H_
