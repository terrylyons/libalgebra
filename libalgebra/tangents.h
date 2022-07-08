//
// Created by user on 07/07/22.
//

#ifndef LIBALGEBRA_LIBALGEBRA_TANGENTS_H_
#define LIBALGEBRA_LIBALGEBRA_TANGENTS_H_

#include <libalgebra/implementation_types.h>
#include <libalgebra/vectors/vectors.h>
#include <libalgebra/algebra.h>
#include <libalgebra/tensor.h>

#include <iosfwd>
#include <type_traits>
#include <utility>


#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
#include <boost/serialization/serialization.hpp>
#endif

namespace alg {


template <typename Vector, typename Fibre = Vector>
class vector_bundle : public Vector
{
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

    vector_bundle() = default;

    explicit vector_bundle(const Vector& point) : Vector(point), m_fibre()
    {}

    explicit vector_bundle(Vector&& point) : Vector(std::move(point)), m_fibre()
    {}

    vector_bundle(Vector&& point, Fibre&& tangent)
        : Vector(std::move(point)), m_fibre(std::move(tangent))
    {}

    vector_bundle(Vector&& point, const Fibre& fibre)
        : Vector(std::move(point)), m_fibre(fibre)
    {}

    vector_bundle(const Vector& point, const Fibre& tangent)
        : Vector(point), m_fibre(tangent)
    {}

    template<typename... Args>
    explicit vector_bundle(Vector&& point, Args&&... args)
        : Vector(std::move(point)), m_fibre(std::forward<Args>(args)...)
    {}

    template<typename... Args>
    explicit vector_bundle(const Vector& point, Args&&... args)
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

    vector_bundle& add_scal_prod(const vector_bundle& rhs, const scalar_type& s);
    vector_bundle& sub_scal_prod(const vector_bundle& rhs, const scalar_type& s);
    vector_bundle& add_scal_div(const vector_bundle& rhs, const rational_type& s);
    vector_bundle& sub_scal_div(const vector_bundle& rhs, const rational_type& s);

    vector_bundle& add_scal_prod(const Vector& rhs, const scalar_type& s);
    vector_bundle& sub_scal_prod(const Vector& rhs, const scalar_type& s);
    vector_bundle& add_scal_div(const Vector& rhs, const rational_type& s);
    vector_bundle& sub_scal_div(const Vector& rhs, const rational_type& s);


//    typename std::enable_if<
//            !std::is_same<
//                    scalar_type,
//                    fibre_scalar_type>::value,
//            vector_bundle&>::type
//    add_scal_prod(const vector_bundle& rhs, const fibre_scalar_type& s);
//
//    typename std::enable_if<
//            !std::is_same<
//                    scalar_type,
//                    fibre_scalar_type>::value,
//            vector_bundle&>::type
//    sub_scal_prod(const vector_bundle& rhs, const fibre_scalar_type& s);
//
//    typename std::enable_if<
//            !std::is_same<
//                    rational_type,
//                    fibre_rational_type>::value,
//            vector_bundle&>::type
//    add_scal_div(const vector_bundle& rhs, const fibre_rational_type& s);
//
//    typename std::enable_if<
//            !std::is_same<
//                    rational_type,
//                    fibre_rational_type>::value,
//            vector_bundle&>::type
//    sub_scal_div(const vector_bundle& rhs, const fibre_rational_type& s);

public:
    vector_bundle& mul_scal_prod(const vector_bundle& rhs, const scalar_type& s, DEG depth);
    vector_bundle& mul_scal_prod(const vector_type& rhs, const scalar_type& s, DEG depth);
    vector_bundle& mul_scal_prod(const vector_bundle& rhs, const scalar_type& s);
    vector_bundle& mul_scal_prod(const vector_type& rhs, const scalar_type& s);

    vector_bundle& mul_scal_div(const vector_bundle& rhs, const rational_type& s, DEG depth);
    vector_bundle& mul_scal_div(const vector_type& rhs, const rational_type& s, DEG depth);
    vector_bundle& mul_scal_div(const vector_bundle& rhs, const rational_type& s);
    vector_bundle& mul_scal_div(const vector_type& rhs, const rational_type& s);

};


template <typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator-(const vector_bundle<Vector, Fibre>& arg)
{
    return {-static_cast<const Vector&>(arg), -arg.fibre()};
}

template <typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator+(const vector_bundle<Vector, Fibre>& lhs, const vector_bundle<Vector, Fibre>& rhs)
{
    return {
        static_cast<const Vector&>(lhs) + static_cast<const Vector&>(rhs),
            lhs.fibre() + rhs.fibre()
    };
}

template <typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator+(const vector_bundle<Vector, Fibre>& lhs, const Vector& rhs)
{
    return {static_cast<const Vector&>(lhs) + rhs, lhs.fibre()};
}

template <typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator+(const Vector& lhs, const vector_bundle<Vector, Fibre>& rhs)
{
    return {lhs + static_cast<const Vector&>(rhs), rhs.fibre()};
}

template <typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator-(const vector_bundle<Vector, Fibre>& lhs, const vector_bundle<Vector, Fibre>& rhs)
{
    return {
        static_cast<const Vector&>(lhs) - static_cast<const Vector&>(rhs),
            lhs.fibre() - rhs.fibre()
    };
}

template <typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator-(const vector_bundle<Vector, Fibre>& lhs, const Vector& rhs)
{
    return {static_cast<const Vector&>(lhs) - rhs, lhs.fibre()};
}

template <typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator-(const Vector& lhs, const vector_bundle<Vector, Fibre>& rhs)
{
    return {lhs - static_cast<const Vector&>(rhs), -rhs.fibre()};
}

/*
 * For scalar multiplication, the vector and fibre types might have different scalar types, so we define
 * scalar multiplication with respect to both and let the compiler figure out which one is more permissive.
 * The same is true for rational division.
 */

template <typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator*(const vector_bundle<Vector, Fibre>& lhs, const typename Vector::SCALAR& rhs)
{
    return {static_cast<const Vector&>(lhs)*rhs, lhs.fibre()*rhs};
}

template <typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator*(const typename Vector::SCALAR& lhs, const vector_bundle<Vector, Fibre>& rhs)
{
    return {lhs*static_cast<const Vector&>(rhs), lhs* rhs.fibre()};
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

template <typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator/(const vector_bundle<Vector, Fibre>& lhs, const typename Vector::RATIONAL& rhs)
{
    return {static_cast<const Vector&>(lhs)/rhs, lhs.fibre()/rhs};
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




template <typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator*(
        const vector_bundle<Vector, Fibre>& left,
        const vector_bundle<Vector, Fibre>& right
        )
{
    const auto& lhs_vec = static_cast<const Vector&>(left);
    const auto& rhs_vec = static_cast<const Vector&>(right);
    return vector_bundle<Vector, Fibre>(
            lhs_vec*rhs_vec,
            lhs_vec* right.fibre() + left.fibre()*rhs_vec
            );
}

template <typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator*(
        const vector_bundle<Vector, Fibre>& left,
        const Vector& rhs_vec
        )
{
    const auto& lhs_vec = static_cast<const Vector&>(left);
    return vector_bundle<Vector, Fibre>(
            lhs_vec*rhs_vec,
            left.fibre()*rhs_vec
            );
}

template <typename Vector, typename Fibre>
vector_bundle<Vector, Fibre> operator*(
        const Vector& lhs_vec,
        const vector_bundle<Vector, Fibre>& right
        )
{
    const auto& rhs_vec = static_cast<const Vector&>(right);
    return vector_bundle<Vector, Fibre>(
            lhs_vec * rhs_vec,
            lhs_vec* right.fibre()
            );
}


/*
 * Now all the in-place operations.
 */

template <typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& operator+=(vector_bundle<Vector, Fibre>& lhs, const vector_bundle<Vector, Fibre>& rhs)
{
    static_cast<Vector&>(lhs) += static_cast<const Vector&>(rhs);
    lhs.fibre() += rhs.fibre();
    return lhs;
}

template <typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& operator+=(vector_bundle<Vector, Fibre>& lhs, const Vector& rhs)
{
    static_cast<Vector&>(lhs) += rhs;
    return lhs;
}

template <typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& operator-=(vector_bundle<Vector, Fibre>& lhs, const vector_bundle<Vector, Fibre>& rhs)
{
    static_cast<Vector&>(lhs) -= static_cast<const Vector&>(rhs);
    lhs.fibre() -= rhs.fibre();
    return lhs;
}

template <typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& operator-=(vector_bundle<Vector, Fibre>& lhs, const Vector& rhs)
{
    static_cast<Vector&>(lhs) -= rhs;
    return lhs;
}

template <typename Vector, typename Fibre>
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

template <typename Vector, typename Fibre>
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


template <typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& operator*=(
        vector_bundle<Vector, Fibre>& lhs,
        const vector_bundle<Vector, Fibre>& rhs
        )
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

template <typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& operator*=(
        vector_bundle<Vector, Fibre>& lhs,
        const Vector& rhs_vec
        )
{
    auto& lhs_vec = static_cast<Vector&>(lhs);

    lhs_vec *= rhs_vec;
    lhs.fibre() *= rhs_vec;
    return lhs;
}




template <typename Vector, typename Fibre>
bool operator==(const vector_bundle<Vector, Fibre>& lhs,
                const vector_bundle<Vector, Fibre>& rhs)
{
    return (static_cast<const Vector&>(lhs) == static_cast<const Vector&>(rhs))
            && (lhs.fibre() == rhs.fibre());
}

template <typename Vector, typename Fibre>
bool operator!=(const vector_bundle<Vector, Fibre>& lhs,
                const vector_bundle<Vector, Fibre>& rhs)
{
    return !operator==(lhs, rhs);
}

template <typename Vector, typename Fibre>
std::ostream& operator<<(std::ostream& os, const vector_bundle<Vector, Fibre>& arg)
{
    return os << '('
              << static_cast<const Vector&>(arg)
              << ", "
              << arg.fibre()
              << ')';
}


template <typename Vector, typename Fibre>
typename std::enable_if<is_algebra<Vector>(), vector_bundle<Vector, Fibre>>::type
commutator(const vector_bundle<Vector, Fibre>& x, const vector_bundle<Vector, Fibre>& y)
{
    vector_bundle<Vector, Fibre> result(x*y);
    result.add_mul(y, x);
    return result;
}




//
//
//template <typename Coeff,
//         DEG Width,
//         DEG Depth,
//         template<typename, typename, typename...> class VectorType,
//         typename... Args>
//class vector_bundle<free_tensor<Coeff, Width, Depth, VectorType, Args...>,
//        free_tensor<Coeff, Width, Depth, VectorType, Args...>>
//    : public free_tensor<Coeff, Width, Depth, VectorType, Args...>
//{
//public:
//
//    using vector_type = free_tensor<Coeff, Width, Depth, VectorType, Args...>;
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
//private:
//
//    fibre_vector_type m_fibre;
//
//public:
//    vector_bundle() = default;
//
//    explicit vector_bundle(const vector_type& point) : vector_type(point), m_fibre()
//    {}
//
//    explicit vector_bundle(vector_type&& point) : vector_type(std::move(point)), m_fibre()
//    {}
//
//    vector_bundle(vector_type&& point, fibre_vector_type&& fibre)
//            : vector_type(std::move(point)), m_fibre(std::move(fibre))
//    {}
//
//    vector_bundle(const vector_type& point, const fibre_vector_type& fibre)
//            : vector_type(point), m_fibre(fibre)
//    {}
//
//    template<typename... CtorArgs>
//    explicit vector_bundle(vector_type&& point, CtorArgs&&... args)
//            : vector_type(std::move(point)), m_fibre(std::forward<Args>(args)...)
//    {}
//
//    template<typename... CtorArgs>
//    explicit vector_bundle(const vector_type& point, CtorArgs&&... args)
//            : vector_type(point), m_fibre(std::forward<CtorArgs>(args)...)
//    {}
//
//    fibre_vector_type& fibre() noexcept
//    {
//        return m_fibre;
//    }
//
//    const fibre_vector_type& fibre() const noexcept
//    {
//        return m_fibre;
//    }
//
//private:
//
//#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
//    friend class boost::serialization::access;
//
//    template<typename Archive>
//    void serialize(Archive& ar, unsigned int const /*version*/)
//    {
//        ar& boost::serialization::base_object<vector_type>(*this);
//        ar& m_fibre;
//    }
//#endif
//
//
//public:
//
////    friend vector_bundle exp(const vector_bundle& arg)
////    {
////        vector_type exponential = exp(static_cast<const vector_type&>(arg));
////        key_type kunit;
////        fibre_vector_type tmp(arg.fibre());
////
////        auto ad_X = [&](const fibre_vector_type& t) {
////            return arg*t - t*arg;
////        };
////
////        scalar_type sign(Depth&1 ? -1 : 1);
////        for (DEG i=Depth; i >= 1; --i) {
////            tmp *= sign / rational_type(i+1);
////            tmp = arg.fibre() + ad_X(tmp);
////            sign *= scalar_type(-1);
////        }
////
////        return vector_bundle(exponential, exponential*tmp);
////    }
//
//    vector_bundle fmexp(const vector_bundle& arg) const
//    {
//        vector_type result(*this), x(arg);
//        key_type kunit;
//
//        const auto& self = static_cast<const vector_type&>(*this);
//
//        auto unit_elt = x.find(kunit);
//        if (unit_elt != x.end() && unit_elt->value() != coeff_type::zero) {
//            x.erase(unit_elt);
//        }
//
//        for (DEG i=Depth; i >= 1; --i) {
//            result.mul_scal_div(x, rational_type(i), Depth-i+1);
//            result += self;
//        }
//
//        return result;
//    }
//
//    vector_bundle& fmexp_inplace(const vector_bundle& arg) const
//    {
//        vector_type self(*this), x(arg);
//
//        key_type kunit;
//
//        auto unit_elt = x.find(kunit);
//        if (unit_elt != x.end() && unit_elt->value() != coeff_type::zero) {
//            x.erase(unit_elt);
//        }
//
//        for (DEG i=Depth; i >= 1; --i) {
//            mul_scal_div(x, rational_type(i), Depth-i+1);
//            *this += self;
//        }
//
//        return *this;
//    }
//
//    vector_bundle& mul_scal_prod(const vector_bundle& rhs, const scalar_type& s, DEG depth);
//    vector_bundle& mul_scal_prod(const vector_type& rhs, const scalar_type& s, DEG depth);
//    vector_bundle& mul_scal_prod(const vector_bundle& rhs, const scalar_type& s);
//    vector_bundle& mul_scal_prod(const vector_type& rhs, const scalar_type& s);
//
//    vector_bundle& mul_scal_div(const vector_bundle& rhs, const rational_type& s, DEG depth);
//    vector_bundle& mul_scal_div(const vector_type& rhs, const rational_type& s, DEG depth);
//    vector_bundle& mul_scal_div(const vector_bundle& rhs, const rational_type& s);
//    vector_bundle& mul_scal_div(const vector_type& rhs, const rational_type& s);
//};


template <typename Coeffs, DEG Width, DEG Depth, template <typename, typename, typename...> class VectorType, typename... Args>
using tensor_bundle = vector_bundle<free_tensor<Coeffs, Width, Depth, VectorType, Args...>>;



/*
 * Below are all the definitions for the extension of vector (algebra) member functions
 * to vector_bundle objects. These are external so they are inherited by default for
 * specializations.
 */







template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& vector_bundle<Vector, Fibre>::add_scal_prod(const vector_bundle& rhs, const scalar_type& s)
{
    vector_type::add_scal_prod(rhs, s);
    m_fibre.add_scal_prod(rhs.fibre(), s);
    return *this;
}
template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& vector_bundle<Vector, Fibre>::sub_scal_prod(const vector_bundle& rhs, const scalar_type& s)
{
    vector_type::sub_scal_prod(rhs, s);
    m_fibre.sub_scal_prod(rhs.fibre(), s);
    return *this;
}

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& vector_bundle<Vector, Fibre>::add_scal_div(const vector_bundle& rhs, const rational_type& s)
{
    vector_type::add_scal_div(rhs, s);
    m_fibre.add_scal_div(rhs.fibre(), s);
    return *this;
}
template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& vector_bundle<Vector, Fibre>::sub_scal_div(const vector_bundle& rhs, const rational_type& s)
{
    vector_type::sub_scal_div(rhs, s);
    m_fibre.sub_scal_div(rhs.fibre(), s);
    return *this;
}

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& vector_bundle<Vector, Fibre>::add_scal_prod(const Vector& rhs, const scalar_type& s)
{
    vector_type::add_scal_prod(rhs, s);
    return *this;
}
template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& vector_bundle<Vector, Fibre>::sub_scal_prod(const Vector& rhs, const scalar_type& s)
{
    vector_type::sub_scal_prod(rhs, s);
    return *this;
}
template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& vector_bundle<Vector, Fibre>::add_scal_div(const Vector& rhs, const rational_type& s)
{
    vector_type::add_scal_div(rhs, s);
    return *this;
}
template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& vector_bundle<Vector, Fibre>::sub_scal_div(const Vector& rhs, const rational_type& s)
{
    vector_type::sub_scal_div(rhs, s);
    return *this;
}

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
template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& vector_bundle<Vector, Fibre>::mul_scal_prod(const vector_type& rhs, const scalar_type& s, DEG depth)
{
    vector_type::mul_scal_prod(rhs, s, depth);
    m_fibre.mul_scal_prod(rhs, s, depth);
    return *this;
}
template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& vector_bundle<Vector, Fibre>::mul_scal_prod(const vector_type& rhs, const scalar_type& s)
{
    vector_type::mul_scal_prod(rhs, s);
    m_fibre.mul_scal_prod(rhs, s);
    return *this;
}
template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& vector_bundle<Vector, Fibre>::mul_scal_div(const vector_type& rhs, const rational_type& s, DEG depth)
{
    vector_type::mul_scal_div(rhs, s, depth);
    m_fibre.mul_scal_div(rhs, s, depth);
    return *this;
}
template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>& vector_bundle<Vector, Fibre>::mul_scal_div(const vector_type& rhs, const rational_type& s)
{
    vector_type::mul_scal_div(rhs, s);
    m_fibre.mul_scal_div(rhs, s);
    return *this;
}

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>&
vector_bundle<Vector, Fibre>::mul_scal_prod(const vector_bundle& rhs, const scalar_type& s, DEG depth)
{
    const auto& rhs_vec = static_cast<const vector_type&>(rhs);
    auto& this_vec = static_cast<vector_type&>(*this);
    vector_type::mul_scal_prod(rhs_vec, s, depth);
    m_fibre = this_vec.mul_scal_prod(rhs.fibre(), s, depth) + fibre().mul_scal_prod(rhs_vec, s, depth);
    return *this;
}
template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>&
vector_bundle<Vector, Fibre>::mul_scal_prod(const vector_bundle& rhs, const scalar_type& s)
{
    const auto& rhs_vec = static_cast<const vector_type&>(rhs);
    auto& this_vec = static_cast<vector_type&>(*this);
    vector_type::mul_scal_prod(rhs_vec, s);
    m_fibre = this_vec.mul_scal_prod(rhs.fibre(), s) + fibre().mul_scal_prod(rhs_vec, s);
    return *this;
}

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>&
vector_bundle<Vector, Fibre>::mul_scal_div(const vector_bundle& rhs, const rational_type& s)
{
    const auto& rhs_vec = static_cast<const vector_type&>(rhs);
    auto& this_vec = static_cast<vector_type&>(*this);
    vector_type::mul_scal_div(rhs_vec, s);
    m_fibre = this_vec.mul_scal_div(rhs.fibre(), s) + fibre().mul_scal_div(rhs_vec, s);
    return *this;
}

template<typename Vector, typename Fibre>
vector_bundle<Vector, Fibre>&
vector_bundle<Vector, Fibre>::mul_scal_div(const vector_bundle& rhs, const rational_type& s, DEG depth)
{
    const auto& rhs_vec = static_cast<const vector_type&>(rhs);
    auto& this_vec = static_cast<vector_type&>(*this);
    vector_type::mul_scal_div(rhs_vec, s, depth);
    m_fibre = this_vec.mul_scal_div(rhs.fibre(), s, depth) + fibre().mul_scal_div(rhs_vec, s, depth);
    return *this;
}

namespace vectors {
namespace dtl {

template <typename Vector, typename Fibre>
struct disable_vector_operations_definition<vector_bundle<Vector, Fibre>>
    : std::true_type
{};

}
}

} // namespace alg

#endif//LIBALGEBRA_LIBALGEBRA_TANGENTS_H_
