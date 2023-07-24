//
// Created by sam on 09/02/2021.
//

#ifndef LIBALGEBRA_BASIS_H
#define LIBALGEBRA_BASIS_H

#include "implementation_types.h"

#include "tags.h"
#include <algorithm>
#include <type_traits>

namespace alg {
namespace basis {

template<typename Basis1, typename Basis2>
struct related_to : std::is_same<Basis1, Basis2> {
};


namespace dtl {

template <typename Basis, typename DegreeTag=typename Basis::degree_tag>
struct basis_trait_impl;

struct resize_info {
    DIMN dimension;
    DIMN size;
    DEG degree;
};

template <typename Basis, DEG D>
struct basis_trait_impl<Basis, with_degree<D>>
{

    template <typename S>
    static void fill_dense_vector(const Basis& basis, const S* begin, const S* end, S* buffer, DEG /*degree*/=0) noexcept
    {
        std::copy(begin, end, buffer);
    }

    static DIMN adjust_dimension(const Basis& basis, DIMN size, DEG degree=0) noexcept
    {
        return basis.start_of_degree(degree);
    }

    static resize_info next_resize_dimension(const Basis& basis, DIMN target, DEG degree=0) noexcept
    {
        if (target == 0 && degree > 0) {
            auto size = basis.start_of_degree(degree+1);
            return {size, size, degree};
        }

        DIMN size;
        while ((size = basis.start_of_degree(degree+1)) < target && degree <= D) {
            ++degree;
        }
        return {size, size, degree};
    }

    static resize_info key_resize_dimension(const Basis& basis, const typename Basis::KEY& key) noexcept
    {
        auto degree = basis.degree(key);
        auto size = basis.start_of_degree(degree+1);
        return {size, size, degree};
    }

    static DIMN size(const Basis& basis, DEG degree=0) noexcept
    { return basis.start_of_degree(degree+1); }

    static DIMN max_dimension(const Basis& basis) noexcept
    { return basis.size(D); }


};

template <typename Basis>
struct basis_trait_impl<Basis, without_degree>
{
    template <typename S>
    static void fill_dense_vector(const Basis& basis, const S* begin, const S* end, S* buffer, DEG /*degree*/=0) noexcept
    {
        std::copy(begin, end, buffer);
    }

    static DIMN size(const Basis& basis, DEG /*degree*/=0) noexcept
    { return basis.size(); }

    static DIMN adjust_dimension(const Basis& basis, DIMN size, DEG /*degree*/=0) noexcept
    { return size; }

    static resize_info next_resize_dimension(const Basis& basis, DIMN target, DEG /*degree*/=0) noexcept
    {
        return {target, target, 0};
    }

    static resize_info key_resize_dimension(const Basis& basis, const typename Basis::KEY& key) noexcept
    {
        auto size = basis.key_to_index(key) + 1;
        return {size, size, 0};
    }

    static DIMN max_dimension(const Basis& basis) noexcept
    { return basis.size(); }

};


} // namespace dtl


template<typename Basis>
struct basis_traits
    : dtl::basis_trait_impl<Basis>
{
    typedef typename Basis::ordering_tag ordering_tag;
    typedef typename Basis::degree_tag degree_tag;


    static DIMN key_to_index(const Basis& basis, const typename Basis::KEY& key) noexcept
    {
        return basis.key_to_index(key);
    }


};

}// namespace basis
}// namespace alg

#endif// LIBALGEBRA_BASIS_H
