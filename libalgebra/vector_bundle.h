//
// Created by user on 07/07/22.
//

#ifndef LIBALGEBRA_LIBALGEBRA_VECTOR_BUNDLE_H_
#define LIBALGEBRA_LIBALGEBRA_VECTOR_BUNDLE_H_

#include "algebra.h"
#include "implementation_types.h"
#include "tensor.h"
#include "vectors.h"

#include <iosfwd>
#include <type_traits>
#include <utility>

#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
#include <boost/serialization/serialization.hpp>
#endif

namespace alg {

namespace dtl {

template<typename Vector, typename Fibre, typename Derived>
class vector_bundle_base : public Vector
{
protected:
    Fibre m_fibre;

public:
    using vector_type = Vector;
    using fibre_vector_type = Fibre;

    using basis_type = typename Vector::BASIS;
    using coeff_type = typename Vector::coefficient_field;
    using key_type = typename Vector::KEY;
    using scalar_type = typename Vector::SCALAR;
    using rational_type = typename Vector::RATIONAL;

    using fibre_basis_type = typename Fibre::BASIS;
    using fibre_coeff_type = typename Fibre::coefficient_field;
    using fibre_key_type = typename Fibre::KEY;
    using fibre_scalar_type = typename Fibre::SCALAR;
    using fibre_rational_type = typename Fibre::RATIONAL;

    // Legacy declarations
    using BASIS = basis_type;
    using SCALAR = scalar_type;
    using RATIONAL = rational_type;
    using KEY = key_type;
    using coefficient_field = coeff_type;

    static_assert(
            std::is_base_of<vectors::vector<basis_type, coeff_type>, Vector>::value,
            "Vector must be a vector type");
    static_assert(
            std::is_base_of<vectors::vector<fibre_basis_type, fibre_coeff_type>, Fibre>::value,
            "Fibre must be a vector type");
    static_assert(
            std::is_same<
                    fibre_scalar_type,
                    decltype(std::declval<scalar_type>() * std::declval<fibre_scalar_type>())>::value
                    && std::is_same<
                            fibre_scalar_type,
                            decltype(std::declval<fibre_scalar_type>() * std::declval<scalar_type>())>::value,
            "fibre scalar type must be multiplicative with vector scalar type");

    vector_bundle_base() = default;

    vector_bundle_base(const vector_bundle_base&) = default;
    vector_bundle_base(vector_bundle_base&&) noexcept = default;

    explicit vector_bundle_base(const Vector& point) : Vector(point), m_fibre()
    {}

    explicit vector_bundle_base(Vector&& point) : Vector(std::move(point)), m_fibre()
    {}

    vector_bundle_base(Vector&& point, Fibre&& tangent)
        : Vector(std::move(point)), m_fibre(std::move(tangent))
    {}

    vector_bundle_base(Vector&& point, const Fibre& fibre)
        : Vector(std::move(point)), m_fibre(fibre)
    {}

    vector_bundle_base(const Vector& point, const Fibre& tangent)
        : Vector(point), m_fibre(tangent)
    {}

    template<typename Key, typename = typename std::enable_if<std::is_same<Key, key_type>::value && std::is_same<Key, fibre_key_type>::value>>
    explicit vector_bundle_base(Key k)
        : Vector(k), m_fibre(k)
    {}

    template<typename... Args>
    explicit vector_bundle_base(Vector&& point, Args&&... args)
        : Vector(std::move(point)), m_fibre(std::forward<Args>(args)...)
    {}

    template<typename... Args>
    explicit vector_bundle_base(const Vector& point, Args&&... args)
        : Vector(point), m_fibre(std::forward<Args>(args)...)
    {}

    Fibre& fibre() noexcept
    {
        return m_fibre;
    }

    const Fibre& fibre() const noexcept
    {
        return m_fibre;
    }

private:
#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive& ar, unsigned int const /*version*/)
    {
        ar& boost::serialization::base_object<Vector>(*this);
        ar& m_fibre;
    }
#endif
public:
    Derived& add_scal_prod(const vector_bundle_base& rhs, const scalar_type& s);
    Derived& sub_scal_prod(const vector_bundle_base& rhs, const scalar_type& s);
    Derived& add_scal_div(const vector_bundle_base& rhs, const rational_type& s);
    Derived& sub_scal_div(const vector_bundle_base& rhs, const rational_type& s);

    Derived& add_scal_prod(const Vector& rhs, const scalar_type& s);
    Derived& sub_scal_prod(const Vector& rhs, const scalar_type& s);
    Derived& add_scal_div(const Vector& rhs, const rational_type& s);
    Derived& sub_scal_div(const Vector& rhs, const rational_type& s);

    template<typename OtherVector>
    typename std::enable_if<!std::is_base_of<vector_bundle_base<Vector, Fibre, Derived>, OtherVector>::value, Derived>::type&
    add_scal_prod(const OtherVector& rhs, const scalar_type& s);

    Derived& mul_scal_prod(const vector_bundle_base& rhs, const scalar_type& s, DEG depth);
    Derived& mul_scal_prod(const vector_type& rhs, const scalar_type& s, DEG depth);
    Derived& mul_scal_prod(const vector_bundle_base& rhs, const scalar_type& s);
    Derived& mul_scal_prod(const vector_type& rhs, const scalar_type& s);

    Derived& mul_scal_div(const vector_bundle_base& rhs, const rational_type& s, DEG depth);
    Derived& mul_scal_div(const vector_type& rhs, const rational_type& s, DEG depth);
    Derived& mul_scal_div(const vector_bundle_base& rhs, const rational_type& s);
    Derived& mul_scal_div(const vector_type& rhs, const rational_type& s);
};

}// namespace dtl

template<typename Vector, typename Fibre = Vector>
class vector_bundle : public dtl::vector_bundle_base<Vector, Fibre, vector_bundle<Vector, Fibre>>
{
public:
    using bundle_base = dtl::vector_bundle_base<Vector, Fibre, vector_bundle>;

    using bundle_base::bundle_base;

private:
#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive& ar, unsigned int const /*version*/)
    {
        ar& boost::serialization::base_object<bundle_base>(*this);
    }
#endif
};

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator-(const vector_bundle<Vector, Fibre>& arg)
{
    return {-static_cast<const Vector&>(arg), -arg.fibre()};
}

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator+(const vector_bundle<Vector, Fibre>& lhs, const vector_bundle<Vector, Fibre>& rhs)
{
    return {
            static_cast<const Vector&>(lhs) + static_cast<const Vector&>(rhs),
            lhs.fibre() + rhs.fibre()};
}

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator+(const vector_bundle<Vector, Fibre>& lhs, const Vector& rhs)
{
    return {static_cast<const Vector&>(lhs) + rhs, lhs.fibre()};
}

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator+(const Vector& lhs, const vector_bundle<Vector, Fibre>& rhs)
{
    return {lhs + static_cast<const Vector&>(rhs), rhs.fibre()};
}

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator-(const vector_bundle<Vector, Fibre>& lhs, const vector_bundle<Vector, Fibre>& rhs)
{
    return {
            static_cast<const Vector&>(lhs) - static_cast<const Vector&>(rhs),
            lhs.fibre() - rhs.fibre()};
}

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator-(const vector_bundle<Vector, Fibre>& lhs, const Vector& rhs)
{
    return {static_cast<const Vector&>(lhs) - rhs, lhs.fibre()};
}

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator-(const Vector& lhs, const vector_bundle<Vector, Fibre>& rhs)
{
    return {lhs - static_cast<const Vector&>(rhs), -rhs.fibre()};
}

/*
 * For scalar multiplication, the vector and fibre types might have different scalar types, so we define
 * scalar multiplication with respect to both and let the compiler figure out which one is more permissive.
 * The same is true for rational division.
 */

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator*(const vector_bundle<Vector, Fibre>& lhs, const typename Vector::SCALAR& rhs)
{
    return {static_cast<const Vector&>(lhs) * rhs, lhs.fibre() * rhs};
}

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator*(const typename Vector::SCALAR& lhs, const vector_bundle<Vector, Fibre>& rhs)
{
    return {lhs * static_cast<const Vector&>(rhs), lhs * rhs.fibre()};
}

//template <typename Vector, typename Fibre>
//typename std::enable_if<
//        !std::is_same<
//                typename Vector::SCALAR,
//                typename Fibre::SCALAR
//        >::value,
//        vector_bundle<Vector, Fibre>
//>::type
//operator*(const vector_bundle<Vector, Fibre>& lhs, const typename Fibre::SCALAR& rhs)
//{
//    return {static_cast<const Vector&>(lhs)*rhs, lhs.fibre()*rhs};
//}

//template <typename Vector, typename Fibre>
//typename std::enable_if <
//        !std::is_same<
//                typename Vector::SCALAR,
//                typename Fibre::SCALAR
//        >::value,
//        vector_bundle<Vector, Fibre>
//>::type
//operator*(const typename Fibre::SCALAR& lhs, const vector_bundle<Vector, Fibre>& rhs)
//{
//    return {lhs*static_cast<const Vector&>(rhs), lhs*rhs.fibre()};
//}

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator/(const vector_bundle<Vector, Fibre>& lhs, const typename Vector::RATIONAL& rhs)
{
    return {static_cast<const Vector&>(lhs) / rhs, lhs.fibre() / rhs};
}

//template <typename Vector, typename Fibre>
//typename std::enable_if<
//        !std::is_same<
//                typename Vector::RATIONAL,
//                typename Fibre::RATIONAL
//        >::value,
//        vector_bundle<Vector, Fibre>
//>::type
//operator/(const vector_bundle<Vector, Fibre>& lhs, const typename Fibre::RATIONAL& rhs)
//{
//    return {static_cast<const Vector&>(lhs)/rhs, lhs.fibre()/rhs};
//}

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator*(
        const vector_bundle<Vector, Fibre>& left,
        const vector_bundle<Vector, Fibre>& right)
{
    const auto& lhs_vec = static_cast<const Vector&>(left);
    const auto& rhs_vec = static_cast<const Vector&>(right);
    return vector_bundle<Vector, Fibre>(
            lhs_vec * rhs_vec,
            lhs_vec * right.fibre() + left.fibre() * rhs_vec);
}

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator*(
        const vector_bundle<Vector, Fibre>& left,
        const Vector& rhs_vec)
{
    const auto& lhs_vec = static_cast<const Vector&>(left);
    return vector_bundle<Vector, Fibre>(
            lhs_vec * rhs_vec,
            left.fibre() * rhs_vec);
}

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator*(
        const Vector& lhs_vec,
        const vector_bundle<Vector, Fibre>& right)
{
    const auto& rhs_vec = static_cast<const Vector&>(right);
    return vector_bundle<Vector, Fibre>(
            lhs_vec * rhs_vec,
            lhs_vec * right.fibre());
}

/*
 * Now all the in-place operations.
 */

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& operator+=(vector_bundle<Vector, Fibre>& lhs, const vector_bundle<Vector, Fibre>& rhs)
{
    static_cast<Vector&>(lhs) += static_cast<const Vector&>(rhs);
    lhs.fibre() += rhs.fibre();
    return lhs;
}

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& operator+=(vector_bundle<Vector, Fibre>& lhs, const Vector& rhs)
{
    static_cast<Vector&>(lhs) += rhs;
    return lhs;
}

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& operator-=(vector_bundle<Vector, Fibre>& lhs, const vector_bundle<Vector, Fibre>& rhs)
{
    static_cast<Vector&>(lhs) -= static_cast<const Vector&>(rhs);
    lhs.fibre() -= rhs.fibre();
    return lhs;
}

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& operator-=(vector_bundle<Vector, Fibre>& lhs, const Vector& rhs)
{
    static_cast<Vector&>(lhs) -= rhs;
    return lhs;
}

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& operator*=(vector_bundle<Vector, Fibre>& lhs, const typename Vector::SCALAR& rhs)
{
    static_cast<Vector&>(lhs) *= rhs;
    lhs.fibre() *= rhs;
    return lhs;
}

//template <typename Vector, typename Fibre>
//typename std::enable_if<
//        !std::is_same<
//                typename Vector::SCALAR,
//                typename Fibre::SCALAR
//        >::value,
//        vector_bundle<Vector, Fibre>&
//>::type
//operator*=(vector_bundle<Vector, Fibre>& lhs, const typename Fibre::SCALAR& rhs)
//{
//    static_cast<Vector&>(lhs) *= rhs;
//    lhs.fibre() *= rhs;
//    return lhs;
//}

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& operator/=(vector_bundle<Vector, Fibre>& lhs, const typename Vector::RATIONAL& rhs)
{
    static_cast<Vector&>(lhs) /= rhs;
    lhs.fibre() /= rhs;
    return lhs;
}

//template <typename Vector, typename Fibre>
//typename std::enable_if<
//        !std::is_same<
//                typename Vector::RATIONAL,
//                typename Fibre::RATIONAL
//        >::value,
//        vector_bundle<Vector, Fibre>&
//>::type
//operator/=(vector_bundle<Vector, Fibre>& lhs, const typename Fibre::RATIONAL& rhs)
//{
//    static_cast<Vector&>(lhs) /= rhs;
//    lhs.fibre() /= rhs;
//    return lhs;
//}

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& operator*=(
        vector_bundle<Vector, Fibre>& lhs,
        const vector_bundle<Vector, Fibre>& rhs)
{
    auto& lhs_vec = static_cast<Vector&>(lhs);
    const auto& rhs_vec = static_cast<const Vector&>(rhs);

    lhs.fibre() *= rhs_vec;
    lhs.fibre() += lhs_vec * rhs.fibre();

    // Do this operation last because otherwise it messes with the
    // calculation of the fibre.
    lhs_vec *= rhs_vec;
    return lhs;
}

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& operator*=(
        vector_bundle<Vector, Fibre>& lhs,
        const Vector& rhs_vec)
{
    auto& lhs_vec = static_cast<Vector&>(lhs);

    lhs_vec *= rhs_vec;
    lhs.fibre() *= rhs_vec;
    return lhs;
}

template<typename Vector, typename Fibre>
bool operator==(const vector_bundle<Vector, Fibre>& lhs,
                const vector_bundle<Vector, Fibre>& rhs)
{
    return (static_cast<const Vector&>(lhs) == static_cast<const Vector&>(rhs))
            && (lhs.fibre() == rhs.fibre());
}

template<typename Vector, typename Fibre>
bool operator!=(const vector_bundle<Vector, Fibre>& lhs,
                const vector_bundle<Vector, Fibre>& rhs)
{
    return !operator==(lhs, rhs);
}

template<typename Vector, typename Fibre>
std::ostream& operator<<(std::ostream& os, const vector_bundle<Vector, Fibre>& arg)
{
    return os << '('
              << static_cast<const Vector&>(arg)
              << ", "
              << arg.fibre()
              << ')';
}

template<typename Vector, typename Fibre>
typename std::enable_if<is_algebra<Vector>(), vector_bundle<Vector, Fibre>>::type
commutator(const vector_bundle<Vector, Fibre>& x, const vector_bundle<Vector, Fibre>& y)
{
    vector_bundle<Vector, Fibre> result(x * y);
    result.add_mul(y, x);
    return result;
}
//
//template<typename Coeff,
//         DEG Width,
//         DEG Depth,
//         template<typename, typename, typename...> class VectorType,
//         typename... Args>
//class vector_bundle<free_tensor<Coeff, Width, Depth, VectorType, Args...>,
//                    free_tensor<Coeff, Width, Depth, VectorType, Args...>>
//    //    : public free_tensor<Coeff, Width, Depth, VectorType, Args...>
//    : public dtl::vector_bundle_base<free_tensor<Coeff, Width, Depth, VectorType, Args...>,
//                                     free_tensor<Coeff, Width, Depth, VectorType, Args...>,
//                                     vector_bundle<free_tensor<Coeff, Width, Depth, VectorType, Args...>,
//                                                   free_tensor<Coeff, Width, Depth, VectorType, Args...>>>
//{
//public:
//    using vector_type = free_tensor<Coeff, Width, Depth, VectorType, Args...>;
//
//    using bundle_base = dtl::vector_bundle_base<vector_type, vector_type, vector_bundle>;
//
//    using fibre_vector_type = vector_type;
//
//    using basis_type = typename vector_type::BASIS;
//    using coeff_type = Coeff;
//    using key_type = typename vector_type::KEY;
//    using scalar_type = typename vector_type::SCALAR;
//    using rational_type = typename vector_type::RATIONAL;
//
//    using fibre_basis_type = basis_type;
//    using fibre_coeff_type = coeff_type;
//    using fibre_key_type = key_type;
//    using fibre_scalar_type = scalar_type;
//    using fibre_rational_type = rational_type;
//
//    // Legacy declarations
//    using BASIS = basis_type;
//    using SCALAR = scalar_type;
//    using RATIONAL = rational_type;
//    using KEY = key_type;
//    using coefficient_field = coeff_type;
//
//    using bundle_base::bundle_base;
//
//#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
//    friend class boost::serialization::access;
//
//    template<typename Archive>
//    void serialize(Archive& ar, unsigned int const /*version*/)
//    {
//        ar& boost::serialization::base_object<bundle_base>(*this);
//    }
//#endif
//
//    vector_bundle fmexp(const vector_bundle& arg) const
//    {
//        vector_bundle result(*this), x(arg);
//        key_type kunit;
//
//        //        const auto& self = static_cast<const vector_type&>(*this);
//
//        auto unit_elt = x.find(kunit);
//        if (unit_elt != x.end() && unit_elt->value() != coeff_type::zero) {
//            x.erase(unit_elt);
//        }
//
//        for (DEG i = Depth; i >= 1; --i) {
//            result.mul_scal_div(x, rational_type(i), Depth - i + 1);
//            result += *this;
//        }
//
//        return result;
//    }
//
//    vector_bundle fmexp(const vector_type& arg) const
//    {
//        vector_bundle result(*this);
//        vector_type x(arg);
//
//        key_type kunit;
//        auto unit_elt = x.find(kunit);
//        if (unit_elt != x.end() && unit_elt->value() != coeff_type::zero) {
//            x.erase(unit_elt);
//        }
//
//        for (DEG i = Depth; i >= 1; --i) {
//            result.mul_scal_div(x, rational_type(i), Depth - i + 1);
//            result += *this;
//        }
//
//        return result;
//    }
//
//    vector_bundle& fmexp_inplace(const vector_bundle& arg)
//    {
//        vector_bundle self(*this), x(arg);
//
//        key_type kunit;
//
//        auto unit_elt = x.find(kunit);
//        if (unit_elt != x.end() && unit_elt->value() != coeff_type::zero) {
//            x.erase(unit_elt);
//        }
//
//        for (DEG i = Depth; i >= 1; --i) {
//            bundle_base::mul_scal_div(x, rational_type(i), Depth - i + 1);
//            *this += self;
//        }
//
//        return *this;
//    }
//
//    vector_bundle& fmexp_inplace(const vector_type& arg)
//    {
//        vector_bundle self(*this);
//        vector_type x(arg);
//
//        key_type kunit;
//        auto unit_elt = x.find(kunit);
//        if (unit_elt != x.end() && unit_elt->value() != coeff_type::zero) {
//            x.erase(unit_elt);
//        }
//
//        for (DEG i = Depth; i >= 1; --i) {
//            bundle_base::mul_scal_div(x, rational_type(i), Depth - i + 1);
//            *this += self;
//        }
//
//        return *this;
//    }
//
//    friend vector_bundle antipode(const vector_bundle& arg)
//    {
//        return {antipode(static_cast<const vector_type&>(arg)), antipode(arg.fibre())};
//    }
//};

template<typename Coeffs, DEG Width, DEG Depth, template<typename, typename, typename...> class VectorType, typename... Args>
using tensor_bundle = vector_bundle<free_tensor<Coeffs, Width, Depth, VectorType, Args...>>;

/*
 * Below are all the definitions for the extension of vector (algebra) member functions
 * to vector_bundle objects. These are external so they are inherited by default for
 * specializations.
 */

//
//
//template<typename Vector, typename Fibre>
//vector_bundle<Vector, Fibre>& vector_bundle<Vector, Fibre>::add_scal_prod(const vector_bundle& rhs, const scalar_type& s)
//{
//    vector_type::add_scal_prod(rhs, s);
//    m_fibre.add_scal_prod(rhs.fibre(), s);
//    return *this;
//}
//template<typename Vector, typename Fibre>
//vector_bundle<Vector, Fibre>& vector_bundle<Vector, Fibre>::sub_scal_prod(const vector_bundle& rhs, const scalar_type& s)
//{
//    vector_type::sub_scal_prod(rhs, s);
//    m_fibre.sub_scal_prod(rhs.fibre(), s);
//    return *this;
//}
//
//template<typename Vector, typename Fibre>
//vector_bundle<Vector, Fibre>& vector_bundle<Vector, Fibre>::add_scal_div(const vector_bundle& rhs, const rational_type& s)
//{
//    vector_type::add_scal_div(rhs, s);
//    m_fibre.add_scal_div(rhs.fibre(), s);
//    return *this;
//}
//template<typename Vector, typename Fibre>
//vector_bundle<Vector, Fibre>& vector_bundle<Vector, Fibre>::sub_scal_div(const vector_bundle& rhs, const rational_type& s)
//{
//    vector_type::sub_scal_div(rhs, s);
//    m_fibre.sub_scal_div(rhs.fibre(), s);
//    return *this;
//}
//
//template<typename Vector, typename Fibre>
//vector_bundle<Vector, Fibre>& vector_bundle<Vector, Fibre>::add_scal_prod(const Vector& rhs, const scalar_type& s)
//{
//    vector_type::add_scal_prod(rhs, s);
//    return *this;
//}
//template<typename Vector, typename Fibre>
//vector_bundle<Vector, Fibre>& vector_bundle<Vector, Fibre>::sub_scal_prod(const Vector& rhs, const scalar_type& s)
//{
//    vector_type::sub_scal_prod(rhs, s);
//    return *this;
//}
//template<typename Vector, typename Fibre>
//vector_bundle<Vector, Fibre>& vector_bundle<Vector, Fibre>::add_scal_div(const Vector& rhs, const rational_type& s)
//{
//    vector_type::add_scal_div(rhs, s);
//    return *this;
//}
//template<typename Vector, typename Fibre>
//vector_bundle<Vector, Fibre>& vector_bundle<Vector, Fibre>::sub_scal_div(const Vector& rhs, const rational_type& s)
//{
//    vector_type::sub_scal_div(rhs, s);
//    return *this;
//}
//
//
//template <typename Vector, typename Fibre>
//template<typename OtherVector>
//vector_bundle<Vector, Fibre>&
//vector_bundle<Vector, Fibre>::add_scal_prod(const OtherVector& rhs, const scalar_type& s)
//{
//    vector_type::add_scal_prod(rhs, s);
//    return *this;
//}

//template<typename Vector, typename Fibre>
//typename std::enable_if<
//        !std::is_same<
//                typename Vector::SCALAR,
//                typename Fibre::SCALAR>::value,
//        vector_bundle<Vector, Fibre>&>::type
//vector_bundle<Vector, Fibre>::add_scal_prod(const vector_bundle& rhs, const fibre_scalar_type& s)
//{
//    vector_type::add_scal_prod(rhs, s);
//    m_fibre.add_scal_prod(rhs.fibre(), s);
//    return *this;
//}
//template<typename Vector, typename Fibre>
//typename std::enable_if<
//        !std::is_same<
//                typename Vector::SCALAR,
//                typename Fibre::SCALAR>::value,
//        vector_bundle<Vector, Fibre>&>::type
//vector_bundle<Vector, Fibre>::sub_scal_prod(const vector_bundle& rhs, const fibre_scalar_type& s)
//{
//    vector_type::sub_scal_prod(rhs, s);
//    m_fibre.sub_scal_prod(rhs.fibre(), s);
//    return *this;
//}
//template<typename Vector, typename Fibre>
//typename std::enable_if<
//        !std::is_same<
//                typename Vector::RATIONAL,
//                typename Fibre::RATIONAL>::value,
//        vector_bundle<Vector, Fibre>&>::type
//vector_bundle<Vector, Fibre>::add_scal_div(const vector_bundle& rhs, const fibre_rational_type& s)
//{
//    vector_type::add_scal_div(rhs, s);
//    m_fibre.add_scal_dv(rhs.fibre(), s);
//    return *this;
//}
//template<typename Vector, typename Fibre>
//typename std::enable_if<
//        !std::is_same<
//                typename Vector::RATIONAL,
//                typename Fibre::RATIONAL>::value,
//        vector_bundle<Vector, Fibre>&>::type
//vector_bundle<Vector, Fibre>::sub_scal_div(const vector_bundle& rhs, const fibre_rational_type& s)
//{
//    vector_type::sub_scal_div(rhs, s);
//    m_fibre.sub_scal_div(rhs.fibre(), s);
//    return *this;
//}
//template<typename Vector, typename Fibre>
//vector_bundle<Vector, Fibre>& vector_bundle<Vector, Fibre>::mul_scal_prod(const vector_type& rhs, const scalar_type& s, DEG depth)
//{
//    vector_type::mul_scal_prod(rhs, s, depth);
//    m_fibre.mul_scal_prod(rhs, s, depth);
//    return *this;
//}
//template<typename Vector, typename Fibre>
//vector_bundle<Vector, Fibre>& vector_bundle<Vector, Fibre>::mul_scal_prod(const vector_type& rhs, const scalar_type& s)
//{
//    vector_type::mul_scal_prod(rhs, s);
//    m_fibre.mul_scal_prod(rhs, s);
//    return *this;
//}
//template<typename Vector, typename Fibre>
//vector_bundle<Vector, Fibre>& vector_bundle<Vector, Fibre>::mul_scal_div(const vector_type& rhs, const rational_type& s, DEG depth)
//{
//    vector_type::mul_scal_div(rhs, s, depth);
//    m_fibre.mul_scal_div(rhs, s, depth);
//    return *this;
//}
//template<typename Vector, typename Fibre>
//vector_bundle<Vector, Fibre>& vector_bundle<Vector, Fibre>::mul_scal_div(const vector_type& rhs, const rational_type& s)
//{
//    vector_type::mul_scal_div(rhs, s);
//    m_fibre.mul_scal_div(rhs, s);
//    return *this;
//}
//
//template<typename Vector, typename Fibre>
//vector_bundle<Vector, Fibre>&
//vector_bundle<Vector, Fibre>::mul_scal_prod(const vector_bundle& rhs, const scalar_type& s, DEG depth)
//{
//    const auto& rhs_vec = static_cast<const vector_type&>(rhs);
//    auto& this_vec = static_cast<vector_type&>(*this);
//    m_fibre.mul_scal_prod(rhs_vec, s, depth);
//    m_fibre += this_vec * (rhs.fibre() * s);
//    vector_type::mul_scal_prod(rhs_vec, s, depth);
//    return *this;
//}
//template<typename Vector, typename Fibre>
//vector_bundle<Vector, Fibre>&
//vector_bundle<Vector, Fibre>::mul_scal_prod(const vector_bundle& rhs, const scalar_type& s)
//{
//    const auto& rhs_vec = static_cast<const vector_type&>(rhs);
//    auto& this_vec = static_cast<vector_type&>(*this);
//    m_fibre.mul_scal_prod(rhs_vec, s);
//    m_fibre += this_vec * (rhs.fibre() * s);
//    vector_type::mul_scal_prod(rhs_vec, s);
//    return *this;
//}
//
//template<typename Vector, typename Fibre>
//vector_bundle<Vector, Fibre>&
//vector_bundle<Vector, Fibre>::mul_scal_div(const vector_bundle& rhs, const rational_type& s)
//{
//    const auto& rhs_vec = static_cast<const vector_type&>(rhs);
//    auto& this_vec = static_cast<vector_type&>(*this);
//    m_fibre.mul_scal_div(rhs_vec, s);
//    m_fibre += this_vec * (rhs.fibre() / s);
//    vector_type::mul_scal_div(rhs_vec, s);
//    return *this;
//}
//
//template<typename Vector, typename Fibre>
//vector_bundle<Vector, Fibre>&
//vector_bundle<Vector, Fibre>::mul_scal_div(const vector_bundle& rhs, const rational_type& s, DEG depth)
//{
//    const auto& rhs_vec = static_cast<const vector_type&>(rhs);
//    auto& this_vec = static_cast<vector_type&>(*this);
//    m_fibre.mul_scal_div(rhs_vec, s, depth);
//    m_fibre += this_vec * (rhs.fibre() / s);
//    vector_type::mul_scal_div(rhs_vec, s, depth);
//    return *this;
//}

namespace dtl {

template<typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::add_scal_prod(const vector_bundle_base& rhs, const scalar_type& s)
{
    vector_type::add_scal_prod(rhs, s);
    m_fibre.add_scal_prod(rhs.m_fibre, s);
    return static_cast<Derived&>(*this);
}
template<typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::sub_scal_prod(const vector_bundle_base& rhs, const scalar_type& s)
{
    vector_type::sub_scal_prod(rhs, s);
    m_fibre.sub_scal_prod(rhs.m_fibre, s);
    return static_cast<Derived&>(*this);
}
template<typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::add_scal_div(const vector_bundle_base& rhs, const rational_type& s)
{
    vector_type::add_scal_div(rhs, s);
    m_fibre.add_scal_div(rhs.m_fibre, s);
    return static_cast<Derived&>(*this);
}
template<typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::sub_scal_div(const vector_bundle_base& rhs, const rational_type& s)
{
    vector_type::sub_scal_div(rhs, s);
    m_fibre.sub_scal_div(rhs.m_fibre, s);
    return static_cast<Derived&>(*this);
}

template<typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::add_scal_prod(const Vector& rhs, const scalar_type& s)
{
    vector_type::add_scal_prod(rhs, s);
    return static_cast<Derived&>(*this);
}
template<typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::sub_scal_prod(const Vector& rhs, const scalar_type& s)
{
    vector_type::sub_scal_prod(rhs, s);
    return static_cast<Derived&>(*this);
}
template<typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::add_scal_div(const Vector& rhs, const rational_type& s)
{
    vector_type::add_scal_div(rhs, s);
    return static_cast<Derived&>(*this);
}
template<typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::sub_scal_div(const Vector& rhs, const rational_type& s)
{
    vector_type::sub_scal_div(rhs, s);
    return static_cast<Derived&>(*this);
}
template<typename Vector, typename Fibre, typename Derived>
template<typename OtherVector>
typename std::enable_if<!std::is_base_of<vector_bundle_base<Vector, Fibre, Derived>, OtherVector>::value, Derived>::type&
vector_bundle_base<Vector, Fibre, Derived>::add_scal_prod(const OtherVector& rhs, const scalar_type& s)
{
    vector_type::add_scal_prod(rhs, s);
    return static_cast<Derived&>(*this);
}
template<typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::mul_scal_prod(const vector_bundle_base& rhs, const scalar_type& s, DEG depth)
{
    m_fibre.mul_scal_prod(rhs, s, depth);
    m_fibre += static_cast<Vector&>(*this) * (rhs.m_fibre * s);
    vector_type::mul_scal_prod(rhs, s, depth);
    return static_cast<Derived&>(*this);
}
template<typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::mul_scal_prod(const vector_type& rhs, const scalar_type& s, DEG depth)
{
    m_fibre.mul_scal_prod(rhs, s, depth);
    vector_type::mul_scal_prod(rhs, s, depth);
    return static_cast<Derived&>(*this);
}
template<typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::mul_scal_prod(const vector_bundle_base& rhs, const scalar_type& s)
{
    m_fibre.mul_scal_prod(rhs, s);
    m_fibre += static_cast<Vector&>(*this) * (rhs.m_fibre * s);
    vector_type::mul_scal_prod(rhs, s);
    return static_cast<Derived&>(*this);
}
template<typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::mul_scal_prod(const vector_type& rhs, const scalar_type& s)
{
    vector_type::mul_scal_prod(rhs, s);
    m_fibre.mul_scal_prod(rhs, s);
    return static_cast<Derived&>(*this);
}
template<typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::mul_scal_div(const vector_bundle_base& rhs, const rational_type& s, DEG depth)
{
    m_fibre.mul_scal_div(rhs, s, depth);
    m_fibre += static_cast<Vector&>(*this) * (rhs.m_fibre / s);
    vector_type::mul_scal_div(rhs, s, depth);
    return static_cast<Derived&>(*this);
}
template<typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::mul_scal_div(const vector_type& rhs, const rational_type& s, DEG depth)
{
    vector_type::mul_scal_div(rhs, s, depth);
    m_fibre.mul_scal_div(rhs, s, depth);
    return static_cast<Derived&>(*this);
}
template<typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::mul_scal_div(const vector_bundle_base& rhs, const rational_type& s)
{
    m_fibre.mul_scal_div(rhs, s);
    m_fibre += static_cast<Vector&>(*this) * (rhs.m_fibre / s);
    vector_type::mul_scal_div(rhs, s);
    return static_cast<Derived&>(*this);
}
template<typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::mul_scal_div(const vector_type& rhs, const rational_type& s)
{
    vector_type::mul_scal_div(rhs, s);
    m_fibre.mul_scal_div(rhs, s);
    return static_cast<Derived&>(*this);
}

}// namespace dtl
//
//namespace vectors {
//namespace dtl {
//
//template<typename Vector, typename Fibre>
//struct disable_vector_operations_definition<vector_bundle<Vector, Fibre>>
//    : std::true_type {};
//
//}// namespace dtl
//}// namespace vectors

}// namespace alg

#endif//LIBALGEBRA_LIBALGEBRA_VECTOR_BUNDLE_H_
