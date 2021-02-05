//
// Created by sam on 01/02/2021.
//

#ifndef LIBALGEBRA_VECTOR_H
#define LIBALGEBRA_VECTOR_H

#include "libalgebra/vectors/sparse_vector.h"


namespace alg {
namespace vectors {

template<typename Basis, typename Field>
struct vector_type_selector {
    typedef sparse_vector<Basis, Field> type;
};


/// Main vector interface
/**
 * Main vector interface for libalgebra.
 *
 *
 * @tparam Basis The basis class for the vector to use
 * @tparam Field The coefficient field to use
 * @tparam VectorImpl The underlying vector class to use. Selected automatically
 * based on the vector_type_selector trait.
 */
template<typename Basis,
        typename Coeffs,
        typename VectorImpl = typename vector_type_selector<Basis, Field>::type
        >
class vector : VectorImpl {
public:

    // Type definitions
    typedef Coeffs FIELD;
    typedef Basis BASIS;
    typedef typename FIELD::S SCALAR;
    typedef typename FIELD::S RATIONAL;
    typedef typename BASIS::KEY KEY;

    // Iterator definitions
    typedef typename VectorImpl::iterator iterator;
    typedef typename VectorImpl::const_iterator const_iterator;

protected:

    // The underlying vector type is accessible from derived classes
    // since we might need to access the class directly in order to
    // optimise some operations.
    typedef VectorImpl UnderlyingVectorType;

public:

    // Pull through function definitions from the underlying vector
    using UnderlyingVectorType::begin;
    using UnderlyingVectorType::end;
    using UnderlyingVectorType::find;
    using UnderlyingVectorType::insert;
    using UnderlyingVectorType::erase;
    using UnderlyingVectorType::operator[];
    using UnderlyingVectorType::clear;
    using UnderlyingVectorType::empty;
    using UnderlyingVectorType::size;

    // Pull the static members from the underlying vector class
    using UnderlyingVectorType::one;
    using UnderlyingVectorType::mone;
    using UnderlyingVectorType::zero;
    using UnderlyingVectorType::basis;

    // norms
    using UnderlyingVectorType::NormL1;

    // Utility
    using UnderlyingVectorType::comp;


protected:

    /// Accessor for underlying vector type for derived classes
    UnderlyingVectorType &underlying_vector() { return *this; }

public:

    // Constructors

    /// Default constructor
    /**
     * Create an instance of an empty vector.
     * This element is neutral with respect to + and -.
     */
    vector(void) : UnderlyingVectorType() {}

    /// Copy constructor
    vector(const UnderlyingVectorType &other) : UnderlyingVectorType(other) {}

    /// Unidimensional constructor.
    /**
    * Create a vector with the value corresponding to key k equal
    * to the given coefficient (default +1).
    */
    explicit vector(const KEY &k, const SCALAR &s = one)
            : UnderlyingVectorType(k, s) {}

    /// Copy from other vector type
    /**
     *
     * @tparam F Other field. Scalar types must be convertible to SCALAR.
     * @tparam V Other underlying vector type.
     */
    template<typename F, typename V>
    explicit vector(const vector<BASIS, F, V> &other) : vector() {
        //TODO: Implement me
    }

    /// Create new vector for result of multiplication
    friend vector create_for_mul(const vector &lhs, const vector &rhs) {
        return UnderlyingVectorType::create_for_mul(lhs, rhs);
    }

public:

    // Swap
    void swap(vector& rhs) {
        UnderlyingVectorType::swap(rhs);
    }


public:

    // Arithmetic
    /*
     * Inplace operations are defined first, then the remainder are derived
     * using the __DECLARE_BINARY_OPERATOR macro.
     */

    /// Additive inverse
    vector operator-(void) const {
        return vector(UnderlyingVectorType::operator-());
    }

    /// Inplace scalar multiply
    vector &operator*=(const SCALAR &s) {
        UnderlyingVectorType::operator*=(s);
        return *this;
    }

    /// Inplace rational divide
    vector &operator/=(const RATIONAL &s) {
        UnderlyingVectorType::operator/=(s);
        return *this;
    }

    /// Inplace addition
    vector &operator+=(const vector &rhs) {
        UnderlyingVectorType::operator+=(rhs);
        return *this;
    }

    /// Inplace subtraction
    vector &operator-=(const vector &rhs) {
        UnderlyingVectorType::operator-=(rhs);
        return *this;
    }

    /// Inplace coordinatewise minimum
    vector &operator&=(const vector &rhs) {
        UnderlyingVectorType::operator&=(rhs);
        return *this;
    }

    /// Inplace coordinatewise maximum
    vector &operator|=(const vector &rhs) {
        UnderlyingVectorType::operator|=(rhs);
        return *this;
    }

    /// Scalar multiply
    __DECLARE_BINARY_OPERATOR(vector, *
    , *=, SCALAR);
    /// Rational divide
    __DECLARE_BINARY_OPERATOR(vector,
    /, /=, RATIONAL);
    /// Addition
    __DECLARE_BINARY_OPERATOR(vector,
    +, +=, vector);
    /// Subtraction
    __DECLARE_BINARY_OPERATOR(vector,
    -, -=, vector);
    /// Coordinatewise minimum
    __DECLARE_BINARY_OPERATOR(vector, &
    , &=, vector);
    /// Coordinatewise maximum
    __DECLARE_BINARY_OPERATOR(vector,
    |, |=, vector);

public:

    //Fused add-scalar-multiply and friends

    /// A version of += fused with scalar multiplication
    /**
     * Version of += fused with scalar multiplication where the
     * right hand side is a unidimensional vector.
     *
     * @param rhs
     * @param s
     * @return
     */
    vector &add_scal_prod(const KEY &rhs, const SCALAR &s) {
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
    vector &add_scal_prod(const vector &rhs, const SCALAR &s) {
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
    vector &sub_scal_prod(const KEY &rhs, const SCALAR &s) {
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
    vector &sub_scal_prod(const vector &rhs, const SCALAR &s) {
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
    vector &add_scal_div(const KEY &rhs, const RATIONAL &s) {
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
    vector &add_scal_div(const vector &rhs, const RATIONAL &s) {
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
    vector &sub_scal_div(const KEY &rhs, const RATIONAL &s) {
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
    vector &sub_scal_div(const vector &rhs, const RATIONAL &s) {
        UnderlyingVectorType::sub_scal_div(rhs, s);
        return *this;
    }

public:

    // Comparison operators

    /// Equality operator
    bool operator==(const vector &rhs) const {
        return UnderlyingVectorType::operator==(rhs);
    }

    /// Non-equality operator
    bool operator!=(const vector &rhs) const {
        return !operator==(rhs);
    }

    /// Lexicographic comparison
    bool operator<(const vector &rhs) const {
        return UnderlyingVectorType::operator<(rhs);
    }


public:

    // Display

    inline friend std::ostream &operator<<(std::ostream &os,
                                           const vector &rhs) {
        return (os << (const UnderlyingVectorType &) rhs);
    }

#if 0  // Not yet implemented
public:

    // Serialization access and methods
    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive &ar, const unsigned int /*version*/) {
        ar & boost::serialization::base_object<UnderlyingVectorType>(*this);
    }
#endif

public:

    // Apply transform methods

    /// Buffered apply transform with separate transforms
    /**
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
    void buffered_apply_binary_transform(
            vector &result,
            const vector &rhs,
            KeyTransform key_transform,
            IndexTransform index_transform,
    ) const {
#if 0
        UnderlyingVectorType::triangular_buffered_apply_transform(
                result, rhs, key_transform, index_transform, max_degree);
#endif
    }

    /// Buffered apply transform with only key transform
    template <typename KeyTransform>
    void buffered_apply_binary_transform(
            vector &result,
            const vector& rhs,
            KeyTransform key_transform
    ) const {
#if 0 // not yet implemented

#endif
    }



};


} // namespace vectors
} // namepace alg
#endif //LIBALGEBRA_VECTOR_H
