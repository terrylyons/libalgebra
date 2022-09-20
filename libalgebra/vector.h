﻿//
// Created by sam on 01/02/2021.
//

#ifndef LIBALGEBRA_VECTOR_H
#define LIBALGEBRA_VECTOR_H

#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
#include <boost/serialization/serialization.hpp>
#endif
#include "sparse_vector.h"

namespace alg {
namespace vectors {

template<typename Basis, typename Field>
struct template_vector_type_selector;

template<typename Basis, typename Coeffs, template<typename, typename, typename...> class VectorImpl = template_vector_type_selector<Basis, Coeffs>::template type, typename... Args>
class vector;

namespace dtl {
class vector_base_access
{
public:
    template<typename B, typename C, template<typename, typename, typename...> class Vector, typename... Args>
    static Vector<B, C, Args...>& convert(Vector<B, C, Args...>& arg)
    {
        return arg;
    }

    template<typename B, typename C, template<typename, typename, typename...> class Vector, typename... Args>
    static const Vector<B, C, Args...>& convert(const Vector<B, C, Args...>& arg)
    {
        return arg;
    }

    template<typename Basis, typename Coeffs, template<typename, typename, typename...> class Vector, typename... Args>
    static Vector<Basis, Coeffs, Args...>& convert(vector<Basis, Coeffs, Vector, Args...>& arg)
    {
        return arg;
    }


//    template<typename B, typename C, template<typename, typename, typename...> class Vector, typename... Args>
//    static const Vector<B, C, Args...>& convert(const Vector<B, C, Args...>& arg)
//    {
//        return arg;
//    }

    template<typename Basis, typename Coeffs, template<typename, typename, typename...> class Vector, typename... Args>
    static const Vector<Basis, Coeffs, Args...>& convert(const vector<Basis, Coeffs, Vector, Args...>& arg)
    {
        return arg;
    }
};

template<typename Vector>
struct data_access : public data_access_base<Vector> {
};

}// namespace dtl

/// Main vector interface
/**
 * @brief Main vector interface for libalgebra.
 *
 * Provides a consistent interface to the underlying storage type. Moreover, it provides some additional
 * cross-vector-type interactions. For example, it provides an interface for adding a dense vector to a sparse
 * vector.
 *
 * All vectors in libalgebra instances of this class template. The VectorImpl parameter determines the actual
 * type of the vector. Most operations are simply passed through to the underlying vector type, but some
 * have additional layers of misdirection to, for example, select the correct implementation to use. This is
 * especially the case for the functions that are used by algebra types for implementing multiplication.
 *
 * @tparam Basis The basis class for the vector to use
 * @tparam Field The coefficient field to use
 * @tparam VectorImpl The underlying vector class to use. Selected automatically
 * based on the vector_type_selector trait.
 */
template<typename Basis, typename Coeffs, template<typename, typename, typename...> class VectorImpl, typename... Args>
class vector : VectorImpl<Basis, Coeffs, Args...>
{
protected:
    // The underlying vector type is accessible from derived classes
    // since we might need to access the class directly in order to
    // optimise some operations.
    typedef VectorImpl<Basis, Coeffs, Args...> UnderlyingVectorType;

    friend class dtl::vector_base_access;

public:
    // Type definitions
    typedef Coeffs coefficient_field;
    typedef Basis BASIS;
    typedef typename coefficient_field::S SCALAR;
    typedef typename coefficient_field::Q RATIONAL;
    typedef typename BASIS::KEY KEY;

    // Iterator definitions
    typedef typename UnderlyingVectorType::iterator iterator;
    typedef typename UnderlyingVectorType::const_iterator const_iterator;

public:
    // Pull through function definitions from the underlying vector
    using UnderlyingVectorType::begin;
    using UnderlyingVectorType::end;
    using UnderlyingVectorType::erase;
    using UnderlyingVectorType::find;
    using UnderlyingVectorType::insert;
    using UnderlyingVectorType::operator[];
    using UnderlyingVectorType::clear;
    using UnderlyingVectorType::empty;
    using UnderlyingVectorType::size;

    // Pull the static members from the underlying vector class
    using UnderlyingVectorType::basis;
    using UnderlyingVectorType::mone;
    using UnderlyingVectorType::one;
    using UnderlyingVectorType::zero;

    // norms
    using UnderlyingVectorType::NormL1;
    using UnderlyingVectorType::NormLInf;

    // Utility
    using UnderlyingVectorType::comp;

    using UnderlyingVectorType::degree_tag;

protected:
    /// Accessor for underlying vector type for derived classes
    UnderlyingVectorType& underlying_vector()
    {
        return *this;
    }

public:
    // Constructors

    /// Default constructor
    /**
     * Create an instance of an empty vector.
     * This element is neutral with respect to + and -.
     */
    vector()
        : UnderlyingVectorType()
    {}

    /// Copy constructor
    vector(const vector& other)
        : UnderlyingVectorType(other)
    {}

    /// Move constructor
    vector(vector&& other) noexcept
        : UnderlyingVectorType(std::move(other))
    {}

    /// Construct from underlying vector type
    explicit vector(const UnderlyingVectorType& other)
        : UnderlyingVectorType(other)
    {}

    /// Unidimensional constructor.
    /**
     * Create a vector with the value corresponding to key k equal
     * to the given coefficient (default +1).
     */
    explicit vector(const KEY& k, const SCALAR& s = one)
        : UnderlyingVectorType(k, s)
    {}

    /// Copy from other vector type
    /**
     *
     * @tparam F Other field. Scalar types must be convertible to SCALAR.
     * @tparam V Other underlying vector type.
     */
    template<typename F, template<typename, typename> class V>
    explicit vector(const vector<BASIS, F, V>& other) : UnderlyingVectorType()
    {
        typename vector<BASIS, F, V>::const_iterator cit;
        for (cit(other.begin()); cit != other.end(); ++cit) {
            operator[](cit->key()) = cit->value();
        }
    }

    /**
     * @brief Construct from iterator
     * @tparam InputIt
     * @param begin
     * @param end
     */
    template<typename InputIt>
    vector(InputIt begin, InputIt end)
        : UnderlyingVectorType()
    {
        UnderlyingVectorType::insert(begin, end);
    }

    /**
     * @brief Construct from pointer to data.
     * @param begin
     * @param end
     */
    vector(SCALAR const* begin, SCALAR const* end)
        : UnderlyingVectorType(begin, end)
    {}

    vector(SCALAR* begin, SCALAR* end) : UnderlyingVectorType(begin, end)
    {}

    vector(DIMN offset, SCALAR const* begin, SCALAR const* end)
        : UnderlyingVectorType(offset, begin, end)
    {}

    vector(DIMN offset, SCALAR* begin, SCALAR* end)
        : UnderlyingVectorType(offset, begin, end)
    {}

    vector& operator=(const vector& other) = default;
    vector& operator=(vector&& other) noexcept = default;

    UnderlyingVectorType& base_vector() noexcept { return *this; }
    const UnderlyingVectorType& base_vector() const noexcept { return *this; }

protected:
    bool ensure_sized_for_degree(const DEG deg)
    {
        return UnderlyingVectorType::ensure_sized_for_degree(deg);
    }

public:
    /// Swap the data in this instance with another
    void swap(vector& rhs)
    {
        UnderlyingVectorType::swap(rhs);
    }

public:
    // Fused add-scalar-multiply and friends

    /// A version of += fused with scalar multiplication
    /**
     * Version of += fused with scalar multiplication where the
     * right hand side is a unidimensional vector.
     *
     * @param rhs
     * @param s
     * @return
     */
    vector& add_scal_prod(const KEY& rhs, const SCALAR& s)
    {
        UnderlyingVectorType::add_scal_prod(rhs, s);
        return *this;
    }

    /// A version of += fused with scalar multiplication
    /**
     * Version of += fused with scalar multiplication where the
     * right hand side is a vector.
     *
     * @param rhs
     * @param s
     * @return
     */
    vector& add_scal_prod(const vector& rhs, const SCALAR& s)
    {
        UnderlyingVectorType::add_scal_prod(rhs, s);
        return *this;
    }

    /// A version of -= fused with scalar multiplication
    /**
     * Version of -= fused with scalar multiplication where the
     * right hand side is a unidimensional vector.
     *
     * @param rhs
     * @param s
     * @return
     */
    vector& sub_scal_prod(const KEY& rhs, const SCALAR& s)
    {
        UnderlyingVectorType::sub_scal_prod(rhs, s);
        return *this;
    }

    /// A version of -= fused with scalar multiplication
    /**
     * Version of -= fused with scalar multiplication where the
     * right hand side is a vector.
     *
     * @param rhs
     * @param s
     * @return
     */
    vector& sub_scal_prod(const vector& rhs, const SCALAR& s)
    {
        UnderlyingVectorType::sub_scal_prod(rhs, s);
        return *this;
    }

    /// A version of += fused with rational division
    /**
     * Version of += fused with rational division where the
     * right hand side is a unidimensional vector.
     *
     * @param rhs
     * @param s
     * @return
     */
    vector& add_scal_div(const KEY& rhs, const RATIONAL& s)
    {
        UnderlyingVectorType::add_scal_div(rhs, s);
        return *this;
    }

    /// A version of += fused with rational division
    /**
     * Version of += fused with rational division where the
     * right hand side is a vector.
     *
     * @param rhs
     * @param s
     * @return
     */
    vector& add_scal_div(const vector& rhs, const RATIONAL& s)
    {
        UnderlyingVectorType::add_scal_div(rhs, s);
        return *this;
    }

    /// A version of -= fused with rational division
    /**
     * Version of -= fused with rational division where the
     * right hand side is a unidimensional vector.
     *
     * @param rhs
     * @param s
     * @return
     */
    vector& sub_scal_div(const KEY& rhs, const RATIONAL& s)
    {
        UnderlyingVectorType::sub_scal_div(rhs, s);
        return *this;
    }

    /// A version of -= fused with rational division
    /**
     * Version of -= fused with rational division where the
     * right hand side is a vector.
     *
     * @param rhs
     * @param s
     * @return
     */
    vector& sub_scal_div(const vector& rhs, const RATIONAL& s)
    {
        UnderlyingVectorType::sub_scal_div(rhs, s);
        return *this;
    }

    // Templated versions of the fused operations. For cross-vector type
    // application.

    template<template<typename, typename, typename...> class V>
    vector& add_scal_prod(const vector<Basis, Coeffs, V>& rhs, const SCALAR& s)
    {
        using VectorType = V<Basis, Coeffs>;
        typename VectorType::const_iterator cit;
        for (cit = rhs.begin(); cit != rhs.end(); ++cit) {
            UnderlyingVectorType::add_scal_prod(cit->key(), cit->value() * s);
        }
        return *this;
    }

    template<template<typename, typename, typename...> class V>
    vector& sub_scal_prod(const vector<Basis, Coeffs, V>& rhs, const SCALAR& s)
    {
        using VectorType = V<Basis, Coeffs>;
        typename VectorType::const_iterator cit;
        for (cit = rhs.begin(); cit != rhs.end(); ++cit) {
            UnderlyingVectorType::sub_scal_prod(cit->key(), cit->value() * s);
        }
        return *this;
    }

    template<template<typename, typename, typename...> class V>
    vector& add_scal_div(const vector<Basis, Coeffs, V>& rhs, const RATIONAL& s)
    {
        using VectorType = V<Basis, Coeffs>;
        typename VectorType::const_iterator cit;
        for (cit = rhs.begin(); cit != rhs.end(); ++cit) {
            UnderlyingVectorType::add_scal_prod(cit->key(), cit->value() / s);
        }
        return *this;
    }

    template<template<typename, typename, typename...> class V>
    vector& sub_scal_div(const vector<Basis, Coeffs, V>& rhs, const RATIONAL& s)
    {
        using VectorType = V<Basis, Coeffs>;
        typename VectorType::const_iterator cit;
        for (cit = rhs.begin(); cit != rhs.end(); ++cit) {
            UnderlyingVectorType::sub_scal_prod(cit->key(), cit->value() / s);
        }
        return *this;
    }

public:
    // Comparison operators

    /// Equality operator
    bool operator==(const vector& rhs) const
    {
        return UnderlyingVectorType::operator==(rhs);
    }

    /// Non-equality operator
    bool operator!=(const vector& rhs) const
    {
        return !operator==(rhs);
    }

    /// Lexicographic comparison
    bool operator<(const vector& rhs) const
    {
        return UnderlyingVectorType::operator<(rhs);
    }

public:
    // Display

    /// Print the vector to an output stream
    inline friend std::ostream& operator<<(std::ostream& os, const vector& rhs)
    {
        return (os << (const UnderlyingVectorType&)rhs);
    }

public:
#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
    // Serialization access and methods
    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive& ar, const unsigned int /*version*/)
    {
        ar& boost::serialization::base_object<UnderlyingVectorType>(*this);
    }
#endif
public:
    // Information methods

    /// Get the maximum degree held by this vector
    DEG degree() const
    {
        return degree_impl(degree_tag);
    }

    /// Test if the maximum degree held by this vector is equal to given value
    bool degree_equals(const DEG degree) const
    {
        return UnderlyingVectorType::degree_equals(degree);
    }

private:
    template<DEG D>
    DEG degree_impl(alg::basis::with_degree<D>) const
    {
        return UnderlyingVectorType::degree();
    }

    DEG degree_impl(alg::basis::without_degree) const
    {
        return 0;
    }

public:
    // Apply transform methods

    /**
     * @brief Buffered apply transform with separate transforms
     *
     * Apply transform to paired vectors using a buffer. Apply using different
     * transforms for by-index or by-key chosen by the under data source.
     *
     * @tparam KeyTransform Key transform type
     * @tparam IndexTransform Index transform type
     * @param result Buffer in which to place the result
     * @param rhs Right hand side buffer
     * @param key_transform Transform to apply by keys (sparse elements)
     * @param index_transform Transform to apply by index (dense elements)
     */
    template<typename KeyTransform, typename IndexTransform>
    void
    buffered_apply_binary_transform(vector& result, const vector& rhs, KeyTransform key_transform,
                                    IndexTransform index_transform) const
    {
        buffered_apply_binary_transform(result, rhs, key_transform, index_transform, UnderlyingVectorType::degree_tag);
    }

    /// Buffered apply transform with only key transform
    template<typename KeyTransform>
    void buffered_apply_binary_transform(vector& result, const vector& rhs, KeyTransform key_transform) const
    {
        buffered_apply_binary_transform(result, rhs, key_transform, UnderlyingVectorType::degree_tag);
    }

    /**
     * @brief Unbuffered apply transform with separate transforms
     *
     * Apply transform to paired vectors using a buffer. Apply using different
     * transforms for by-index or by-key chosen by the under data source.
     *
     * @tparam KeyTransform Key transform type
     * @tparam IndexTransform Index transform type
     * @param rhs Right hand side buffer
     * @param key_transform Transform to apply by keys (sparse elements)
     * @param index_transform Transform to apply by index (dense elements)
     */
    template<typename KeyTransform, typename IndexTransform>
    void
    unbuffered_apply_binary_transform(const vector& rhs, KeyTransform key_transform, IndexTransform index_transform)
    {
        unbuffered_apply_binary_transform(rhs, key_transform, index_transform, UnderlyingVectorType::degree_tag);
    }

    /// Buffered apply transform with only key transform
    template<typename KeyTransform>
    void unbuffered_apply_binary_transform(const vector& rhs, KeyTransform key_transform)
    {
        unbuffered_apply_binary_transform(rhs, key_transform, UnderlyingVectorType::degree_tag);
    }

    /**
     * @brief Buffered apply transform with separate transforms up to max degree
     *
     * @tparam KeyTransform Key transform type
     * @tparam IndexTransform Index transform type
     * @param result Buffer in which to place the result
     * @param rhs Right hand side buffer
     * @param key_transform Transform to apply by keys (sparse elements)
     * @param index_transform Transform to apply by index (dense elements)
     * @param max_depth Maximum depth to compute the result
     */
    template<typename KeyTransform, typename IndexTransform>
    void
    buffered_apply_binary_transform(vector& result, const vector& rhs, KeyTransform key_transform,
                                    IndexTransform index_transform, const DEG max_depth) const
    {
        UnderlyingVectorType::triangular_buffered_apply_transform(result, rhs, key_transform, index_transform,
                                                                  max_depth);
    }

    /**
     * @brief Unbuffered apply transform with separate transforms up to max degree
     *
     * @tparam KeyTransform Key transform type
     * @tparam IndexTransform Index transform type
     * @param rhs Right hand side buffer
     * @param key_transform Transform to apply by keys (sparse elements)
     * @param index_transform Transform to apply by index (dense elements)
     * @param max_depth Maximum depth to compute the result
     */
    template<typename KeyTransform, typename IndexTransform>
    void
    unbuffered_apply_binary_transform(const vector& rhs, KeyTransform key_transform, IndexTransform index_transform,
                                      const DEG max_depth)
    {
        UnderlyingVectorType::triangular_unbuffered_apply_binary_transform(rhs, key_transform, index_transform,
                                                                           max_depth);
    }

    /// Buffered apply transform with only key transform up to max degree
    template<typename KeyTransform>
    void
    buffered_apply_binary_transform(vector& result, const vector& rhs, KeyTransform key_transform,
                                    const DEG max_depth) const
    {
        UnderlyingVectorType::triangular_buffered_apply_binary_transform(result, rhs, key_transform, max_depth);
    }


private:
    template<typename KeyTransform>
    void
    buffered_apply_binary_transform(vector& result, const vector& rhs, KeyTransform key_transform,
                                    alg::basis::without_degree) const
    {
        UnderlyingVectorType::square_buffered_apply_binary_transform(result, rhs, key_transform);
    }

    template<DEG D, typename KeyTransform>
    void
    buffered_apply_binary_transform(vector& result, const vector& rhs, KeyTransform key_transform,
                                    alg::basis::with_degree<D>) const
    {
        UnderlyingVectorType::triangular_buffered_apply_binary_transform(result, rhs, key_transform, D);
    }

    template<typename KeyTransform, typename IndexTransform>
    void
    buffered_apply_binary_transform(vector& result, const vector& rhs, KeyTransform key_transform,
                                    IndexTransform index_transform, alg::basis::without_degree) const
    {
        UnderlyingVectorType::square_buffered_apply_binary_transform(result, rhs, key_transform, index_transform);
    }

    template<DEG D, typename KeyTransform, typename IndexTransform>
    void
    buffered_apply_binary_transform(vector& result, const vector& rhs, KeyTransform key_transform,
                                    IndexTransform index_transform, alg::basis::with_degree<D>) const
    {
        UnderlyingVectorType::triangular_buffered_apply_binary_transform(result, rhs, key_transform, index_transform,
                                                                         D);
    }

    template<typename KeyTransform>
    void unbuffered_apply_binary_transform(const vector& rhs, KeyTransform key_transform, alg::basis::without_degree)
    {
        vector result;
        UnderlyingVectorType::square_buffered_apply_binary_transform(result, rhs, key_transform);
        swap(result);
    }

    template<DEG D, typename KeyTransform>
    void unbuffered_apply_binary_transform(const vector& rhs, KeyTransform key_transform, alg::basis::with_degree<D>)
    {
        vector result;
        UnderlyingVectorType::triangular_buffered_apply_binary_transform(result, rhs, key_transform, D);
        swap(result);
    }

    template<typename KeyTransform, typename IndexTransform>
    void
    unbuffered_apply_binary_transform(const vector& rhs, KeyTransform key_transform, IndexTransform index_transform,
                                      alg::basis::without_degree)
    {
        vector result;
        UnderlyingVectorType::square_buffered_apply_binary_transform(result, rhs, key_transform, index_transform);
        swap(result);
    }

    template<DEG D, typename KeyTransform, typename IndexTransform>
    void
    unbuffered_apply_binary_transform(const vector& rhs, KeyTransform key_transform, IndexTransform index_transform,
                                      alg::basis::with_degree<D>)
    {
        UnderlyingVectorType::triangular_unbuffered_apply_binary_transform(rhs, key_transform, index_transform, D);
    }

public:
    // Methods for operator implementation

    /// Apply a transform inplace to the vector with buffering
    template<typename Transform>
    void buffered_apply_unary_transform(vector& result, Transform transform) const
    {
        buffered_apply_unary_transform_impl(result, transform, UnderlyingVectorType::degree_tag);
    }

private:
    template<DEG D, typename Transform>
    void buffered_apply_unary_transform_impl(vector& result, Transform transform, alg::basis::with_degree<D>) const
    {
        UnderlyingVectorType::buffered_apply_unary_transform(result, transform, D);
    }

    template<typename Transform>
    void buffered_apply_unary_transform_impl(vector& result, Transform transform, alg::basis::without_degree) const
    {
        UnderlyingVectorType::buffered_apply_unary_transform(result, transform);
    }


public:

    template <typename Vector>
    friend typename
    std::enable_if<std::is_base_of<vector, Vector>::value, Vector>::type
    operator-(const Vector& arg)
    {
        Vector result;
        UnderlyingVectorType::apply_unary_operation(result, arg, Coeffs::uminus);
        return result;
    }

    template <typename Vector, typename Scalar>
    friend
    typename std::enable_if<
            std::is_base_of<vector, Vector>::value &&
                    std::is_constructible<SCALAR, const Scalar&>::value,
            Vector>::type
    operator*(const Vector& arg, const Scalar& scal)
    {
        Vector result;
        SCALAR s(scal);
        UnderlyingVectorType::apply_unary_operation(result, arg,
                                                    [=](const SCALAR& a) { return Coeffs::mul(a, s); }
        );
        return result;
    }

    template <typename Vector>
    friend
    typename std::enable_if<
            std::is_base_of<vector, Vector>::value,
            Vector>::type
    operator*(const Vector& arg, const SCALAR& scal)
    {
        Vector result;
        UnderlyingVectorType::apply_unary_operation(result, arg,
                                                    [=](const SCALAR& a) { return Coeffs::mul(a, scal); }
        );
        return result;
    }

    template <typename Vector, typename Scalar>
    friend
    typename std::enable_if<
            std::is_base_of<vector, Vector>::value &&
            std::is_constructible<SCALAR, const Scalar&>::value,
            Vector>::type
    operator*(const Scalar& scal, const Vector& arg)
    {
        Vector result;
        SCALAR s(scal);
        UnderlyingVectorType::apply_unary_operation(result, arg,
                    [=](const SCALAR& a) { return Coeffs::mul(s, a); }
                );
        return result;
    }

    template <typename Vector>
    friend
    typename std::enable_if<
            std::is_base_of<vector, Vector>::value,
            Vector>::type
    operator*(const SCALAR& scal, const Vector& arg)
    {
        Vector result;
        UnderlyingVectorType::apply_unary_operation(result, arg,
                    [=](const SCALAR& a) { return Coeffs::mul(scal, a); }
                );
        return result;
    }

    template<typename Vector, typename Rational>
    friend typename std::enable_if<
            std::is_base_of<vector, Vector>::value && std::is_constructible<RATIONAL, const Rational&>::value,
            Vector>::type
    operator/(const Vector& arg, const Rational& scal)
    {
        Vector result;
        RATIONAL r(scal);
        UnderlyingVectorType::apply_unary_operation(result, arg,
                                                    [=](const SCALAR& a) { return Coeffs::div(a, r); });
        return result;
    }

    template<typename Vector1, typename Vector2>
    friend
    typename std::enable_if<
            std::is_base_of<vector, Vector1>::value && std::is_base_of<vector, Vector2>::value,
            Vector1>::type
    operator+(const Vector1& lhs, const Vector2& rhs)
    {
        Vector1 result;
        UnderlyingVectorType::apply_flat_binary_operation(
                result, lhs, rhs, Coeffs::template add<>);
        return result;
    }

    template <typename Vector,
             template <typename, typename, typename...> class OtherVType,
             typename... OtherArgs>
    friend typename std::enable_if<
            std::is_base_of<vector, Vector>::value,
            Vector>::type
    operator+(const Vector& lhs, const vector<Basis, Coeffs, OtherVType, OtherArgs...>& rhs)
    {
        Vector result(lhs);
        for (auto item : rhs) {
            result.add_scal_prod(item.key(), item.value());
        }
        return result;
    }

    template<typename Vector1, typename Vector2>
    friend
    typename std::enable_if<
            std::is_base_of<vector, Vector1>::value && std::is_base_of<vector, Vector2>::value,
            Vector1>::type
    operator-(const Vector1& lhs, const Vector2& rhs)
    {
        Vector1 result;
        UnderlyingVectorType::apply_flat_binary_operation(
                result, lhs, rhs, Coeffs::template sub<>);
        return result;
    }

    template <typename Vector,
             template <typename, typename, typename...> class OtherVType,
             typename... OtherArgs>
    friend typename std::enable_if<
            std::is_base_of<vector, Vector>::value,
            Vector>::type
    operator-(const Vector& lhs, const vector<Basis, Coeffs, OtherVType, OtherArgs...>& rhs)
    {
        Vector result(lhs);
        for (auto item : rhs) {
            result.sub_scal_prod(item.key(), item.value());
        }
        return result;
    }

    template<typename Vector, typename Scalar>
    friend typename std::enable_if<
            std::is_base_of<vector, Vector>::value && std::is_constructible<SCALAR, const Scalar&>::value,
            Vector>::type&
    operator*=(Vector& arg, const Scalar& scal)
    {
        SCALAR s(scal);
        UnderlyingVectorType::apply_inplace_unary_op(arg,
                                                     [=](const SCALAR& a) { return Coeffs::mul(a, s); });
        return arg;
    }

    template<typename Vector, typename Rational>
    friend typename std::enable_if<
            std::is_base_of<vector, Vector>::value && std::is_constructible<RATIONAL, const Rational&>::value,
            Vector>::type&
    operator/=(Vector& arg, const Rational& rat)
    {
        RATIONAL r(rat);
        UnderlyingVectorType::apply_inplace_unary_op(arg,
                                                     [=](const SCALAR& a) { return Coeffs::div(a, r); });
        return arg;
    }

    template <typename Vector1, typename Vector2>
    friend typename std::enable_if<
            std::is_base_of<vector, Vector1>::value &&
            std::is_base_of<vector, Vector2>::value,
            Vector1>::type&
    operator+=(Vector1& lhs, const Vector2& rhs)
    {
        UnderlyingVectorType::apply_inplace_flat_binary_op(lhs, rhs, Coeffs::template add<>);
        return lhs;
    }

    template <typename Vector, template <typename, typename, typename...> class OtherVType, typename... OtherArgs>
    friend typename std::enable_if<std::is_base_of<vector, Vector>::value, Vector>::type&
    operator+=(Vector& lhs, const vector<Basis, Coeffs, OtherVType, OtherArgs...>& rhs)
    {
        for (auto item : rhs) {
            lhs.add_scal_prod(item.key(), item.value());
        }
        return lhs;
    }

    template <typename Vector1, typename Vector2>
    friend typename std::enable_if<
            std::is_base_of<vector, Vector1>::value &&
            std::is_base_of<vector, Vector2>::value,
            Vector1>::type&
    operator-=(Vector1& lhs, const Vector2& rhs)
    {
        UnderlyingVectorType::apply_inplace_flat_binary_op(lhs, rhs, Coeffs::template sub<>);
        return lhs;
    }

    template <typename Vector, template<typename, typename, typename...> class OtherVType, typename... OtherArgs>
    friend typename std::enable_if<std::is_base_of<vector, Vector>::value, Vector>::type&
    operator-=(Vector& lhs, const vector<Basis, Coeffs, OtherVType, OtherArgs...>& rhs)
    {
        for (auto item : rhs) {
            lhs.sub_scal_prod(item.key(), item.value());
        }
        return lhs;
    }


    /*
     * Ordering operators are interesting because not all coefficient rings
     * are necessarily ordered. To make sure things are only defined when
     * appropriate, we can include the order as a SFINAE template parameter
     * to the order operations. For the time being, default to std::less,
     * until proper support is added to the coefficient trait.
     */

    template <typename Vector1, typename Vector2, typename Order=std::less<SCALAR>>
    friend typename std::enable_if<
            std::is_base_of<vector, Vector1>::value &&
            std::is_base_of<vector, Vector2>::value,
            Vector1>::type
    operator|(const Vector1& lhs, const Vector2& rhs)
    {
        Vector1 result;
        Order cmp;
        UnderlyingVectorType::apply_flat_binary_operation(
                result, lhs, rhs, [=](const SCALAR& l, const SCALAR& r)
                { return std::max(l, r, cmp); });
        return result;
    }

    template <typename Vector1, typename Vector2, typename Order=std::less<SCALAR>>
    friend typename std::enable_if<
            std::is_base_of<vector, Vector1>::value &&
            std::is_base_of<vector, Vector2>::value,
            Vector1>::type
    operator&(const Vector1& lhs, const Vector2& rhs)
    {
        Vector1 result;
        Order cmp;
        UnderlyingVectorType::apply_flat_binary_operation(result, lhs, rhs,
                  [=](const SCALAR& l, const SCALAR& r) { return std::min(l, r, cmp); });
        return result;
    }

    template <typename Vector1, typename Vector2, typename Order=std::less<SCALAR>>
    friend typename std::enable_if<
            std::is_base_of<vector, Vector1>::value && std::is_base_of<vector, Vector2>::value,
            Vector1>::type
    operator|=(Vector1& lhs, const Vector2& rhs)
    {
        Order cmp;
        UnderlyingVectorType::apply_inplace_flat_binary_op(lhs, rhs,
                           [=](const SCALAR& l, const SCALAR& r) { return std::max(l ,r, cmp); });
        return lhs;
    }

    template <typename Vector1, typename Vector2, typename Order=std::less<SCALAR>>
    friend typename std::enable_if<
            std::is_base_of<vector, Vector1>::value && std::is_base_of<vector, Vector2>::value,
            Vector1>::type
    operator&=(Vector1& lhs, const Vector2& rhs)
    {
        Order cmp;
        UnderlyingVectorType::apply_inplace_flat_binary_op(lhs, rhs,
                           [=](const SCALAR& l, const SCALAR& r) { return std::min(l ,r, cmp); });
        return lhs;
    }
};

}// namespace vectors

namespace utils {


template<typename T>
struct is_vector_type {
private:

    template <typename B, typename C, template <typename, typename, typename...> class Type, typename... Args>
    static std::true_type check(vectors::vector<B, C, Type, Args...>&);

    static std::false_type check(...);

public:

    static constexpr bool value = decltype(check(std::declval<T>()))::value;

};


} // namespace utils


}// namespace alg
#endif// LIBALGEBRA_VECTOR_H