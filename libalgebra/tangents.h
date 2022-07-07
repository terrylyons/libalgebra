//
// Created by user on 07/07/22.
//

#ifndef LIBALGEBRA_LIBALGEBRA_TANGENTS_H_
#define LIBALGEBRA_LIBALGEBRA_TANGENTS_H_

#include <libalgebra/implementation_types.h>
#include <libalgebra/vectors/vectors.h>
#include <libalgebra/algebra.h>
#include <libalgebra/tensor.h>

#include <type_traits>
#include <utility>

#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
#include <boost/serialization/serialization.hpp>
#endif

namespace alg {


template <typename Vector, typename TangentVector=Vector>
class tangent_vector : public Vector
{
    TangentVector m_tangent;

public:
    using vector_type = Vector;
    using tangent_vector_type = TangentVector;

    using basis_type = typename Vector::BASIS;
    using coeff_type = typename Vector::coefficient_field;
    using key_type = typename Vector::KEY;
    using scalar_type = typename Vector::SCALAR;
    using rational_type = typename Vector::RATIONAL;

    using tangent_basis_type = typename TangentVector::BASIS;
    using tangent_coeff_type = typename TangentVector::coefficient_field;
    using tangent_key_type = typename TangentVector::KEY;
    using tangent_scalar_type = typename TangentVector::SCALAR;
    using tangent_rational_type = typename TangentVector::RATIONAL;

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
            std::is_base_of<vectors::vector<tangent_basis_type, tangent_coeff_type>, TangentVector>::value,
            "TangentVector must be a vector type");
    static_assert(
            std::is_same<
                    tangent_scalar_type,
                    decltype(std::declval<scalar_type>() * std::declval<tangent_scalar_type>())>::value_type
                    && std::is_same<
                            tangent_scalar_type,
                            decltype(std::declval<tangent_scalar_type>() * std::declval<scalar_type>())>::value_type,
            "tangent scalar type must be multiplicative with vector scalar type");

    explicit tangent_vector(const Vector& point) : Vector(point), m_tangent()
    {}

    explicit tangent_vector(Vector&& point) : Vector(std::move(point)), m_tangent()
    {}

    tangent_vector(Vector&& point, TangentVector&& tangent)
        : Vector(std::move(point)), m_tangent(std::move(tangent))
    {}

    tangent_vector(const Vector& point, const TangentVector& tangent)
        : Vector(point), m_tangent(tangent)
    {}

    template<typename... Args>
    explicit tangent_vector(Vector&& point, Args&&... args)
        : Vector(std::move(point)), m_tangent(std::forward<Args>(args)...)
    {}

    template<typename... Args>
    explicit tangent_vector(const Vector& point, Args&&... args)
        : Vector(point), m_tangent(std::forward<Args>(args)...)
    {}

    TangentVector& tangent() noexcept
    {
        return m_tangent;
    }

    const TangentVector& tangent() const noexcept
    {
        return m_tangent;
    }

private:
#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive& ar, unsigned int const /*version*/)
    {
        ar& boost::serialization::base_object<Vector>(*this);
        ar& m_tangent;
    }
#endif

public:
    tangent_vector& mul_scal_prod(const tangent_vector& rhs, const scalar_type& s, DEG depth);
    tangent_vector& mul_scal_prod(const vector_type& rhs, const scalar_type& s, DEG depth);
    tangent_vector& mul_scal_prod(const tangent_vector& rhs, const scalar_type& s);
    tangent_vector& mul_scal_prod(const vector_type& rhs, const scalar_type& s);

    tangent_vector& mul_scal_div(const tangent_vector& rhs, const rational_type& s, DEG depth);
    tangent_vector& mul_scal_div(const vector_type& rhs, const rational_type& s, DEG depth);
    tangent_vector& mul_scal_div(const tangent_vector& rhs, const rational_type& s);
    tangent_vector& mul_scal_div(const vector_type& rhs, const rational_type& s);

};

template <typename Vector, typename TangentVector>
tangent_vector<Vector, TangentVector> operator*(
        const tangent_vector<Vector, TangentVector>& left,
        const tangent_vector<Vector, TangentVector>& right
        )
{
    const auto& lhs_vec = static_cast<const Vector&>(left);
    const auto& rhs_vec = static_cast<const Vector&>(right);
    return tangent_vector<Vector, TangentVector>(
            lhs_vec*rhs_vec,
            lhs_vec*right.tangent() + left.tanget()*rhs_vec
            );
}

template <typename Vector, typename TangentVector>
tangent_vector<Vector, TangentVector> operator*(
        const tangent_vector<Vector, TangentVector>& left,
        const Vector& rhs_vec
        )
{
    const auto& lhs_vec = static_cast<const Vector&>(left);
    return tangent_vector<Vector, TangentVector>(
            lhs_vec*rhs_vec,
            left.tangent()*rhs_vec
            );
}

template <typename Vector, typename TangentVector>
tangent_vector<Vector, TangentVector> operator*(
        const Vector& lhs_vec,
        const tangent_vector<Vector, TangentVector>& right
        )
{
    const auto& rhs_vec = static_cast<const Vector&>(right);
    return tangent_vector<Vector, TangentVector>(
            lhs_vec * rhs_vec,
            lhs_vec*right.tangent()
            );
}

template <typename Vector, typename TangentVector>
tangent_vector<Vector, TangentVector>& operator*=(
        tangent_vector<Vector, TangentVector>& lhs,
        const tangent_vector<Vector, TangentVector>& rhs
        )
{
    const auto& lhs_vec = static_cast<const Vector&>(lhs);
    const auto& rhs_vec = static_cast<const Vector&>(rhs);

    lhs_vec *= rhs_vec;
    (lhs.tangent() *= rhs_vec) += (lhs_vec*rhs.tangent());
    return lhs;
}

template <typename Vector, typename TangentVector>
tangent_vector<Vector, TangentVector>& operator*=(
        tangent_vector<Vector, TangentVector>& lhs,
        const Vector& rhs_vec
        )
{
    const auto& lhs_vec = static_cast<const Vector&>(lhs);

    lhs_vec *= rhs_vec;
    lhs.tangent() *= rhs_vec;
    return lhs;
}








template <typename Vector, typename TangentVector>
typename std::enable_if<is_algebra<Vector>(), tangent_vector<Vector, TangentVector>>::type
commutator(const tangent_vector<Vector, TangentVector>& x, const tangent_vector<Vector, TangentVector>& y)
{
    tangent_vector<Vector, TangentVector> result(x*y);
    result.add_mul(y, x);
    return result;
}






template <typename Coeff,
         DEG Width,
         DEG Depth,
         template<typename, typename, typename...> class VectorType,
         typename... Args>
class tangent_vector<free_tensor<Coeff, Width, Depth, VectorType, Args...>,
        free_tensor<Coeff, Width, Depth, VectorType, Args...>>
    : public free_tensor<Coeff, Width, Depth, VectorType, Args...>
{
public:

    using vector_type = free_tensor<Coeff, Width, Depth, VectorType, Args...>;
    using tangent_vector_type = vector_type;

    using basis_type = typename vector_type::BASIS;
    using coeff_type = Coeff;
    using key_type = typename vector_type::KEY;
    using scalar_type = typename vector_type::SCALAR;
    using rational_type = typename vector_type::RATIONAL;

    using tangent_basis_type = basis_type;
    using tangent_coeff_type = coeff_type;
    using tangent_key_type = key_type;
    using tangent_scalar_type = scalar_type;
    using tangent_rational_type = rational_type;

    // Legacy declarations
    using BASIS = basis_type;
    using SCALAR = scalar_type;
    using RATIONAL = rational_type;
    using KEY = key_type;
    using coefficient_field = coeff_type;

private:

    tangent_vector_type m_tangent;

public:
    explicit tangent_vector(const vector_type& point) : vector_type(point), m_tangent()
    {}

    explicit tangent_vector(vector_type&& point) : vector_type(std::move(point)), m_tangent()
    {}

    tangent_vector(vector_type&& point, tangent_vector_type&& tangent)
            : vector_type(std::move(point)), m_tangent(std::move(tangent))
    {}

    tangent_vector(const vector_type& point, const tangent_vector_type& tangent)
            : vector_type(point), m_tangent(tangent)
    {}

    template<typename... CtorArgs>
    explicit tangent_vector(vector_type&& point, CtorArgs&&... args)
            : vector_type(std::move(point)), m_tangent(std::forward<Args>(args)...)
    {}

    template<typename... CtorArgs>
    explicit tangent_vector(const vector_type& point, CtorArgs&&... args)
            : vector_type(point), m_tangent(std::forward<CtorArgs>(args)...)
    {}

    tangent_vector_type& tangent() noexcept
    {
        return m_tangent;
    }

    const tangent_vector_type& tangent() const noexcept
    {
        return m_tangent;
    }

private:

#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive& ar, unsigned int const /*version*/)
    {
        ar& boost::serialization::base_object<vector_type>(*this);
        ar& m_tangent;
    }
#endif


public:

    friend tangent_vector exp(const tangent_vector& arg)
    {
        vector_type exponential = exp(static_cast<const vector_type&>(arg));
        key_type kunit;
        tangent_vector_type tmp(arg.tangent());

        auto ad_X = [&](const tangent_vector_type& t) {
            return arg*t + t*arg;
        };

        scalar_type sign(Depth&1 ? -1 : 1);
        for (DEG i=Depth; i >= 1; --i) {
            tmp *= sign / rational_type(i+1);
            tmp = arg.tangent() + ad_X(tmp);
        }

        return tangent_vector(exponential, exponential*tmp);
    }

    tangent_vector fmexp(const tangent_vector& arg) const
    {
        vector_type result(*this), x(arg);
        key_type kunit;

        const auto& self = static_cast<const vector_type&>(*this);

        auto unit_elt = x.find(kunit);
        if (unit_elt != x.end() && unit_elt->value() != coeff_type::zero) {
            x.erase(unit_elt);
        }

        for (DEG i=Depth; i >= 1; --i) {
            result.mul_scal_div(x, rational_type(i), Depth-i+1);
            result += self;
        }

        return result;
    }

    tangent_vector& fmexp_inplace(const tangent_vector& arg) const
    {
        vector_type self(*this), x(arg);

        key_type kunit;

        auto unit_elt = x.find(kunit);
        if (unit_elt != x.end() && unit_elt->value() != coeff_type::zero) {
            x.erase(unit_elt);
        }

        for (DEG i=Depth; i >= 1; --i) {
            mul_scal_div(x, rational_type(i), Depth-i+1);
            *this += self;
        }

        return *this;
    }

    tangent_vector& mul_scal_prod(const tangent_vector& rhs, const scalar_type& s, DEG depth);
    tangent_vector& mul_scal_prod(const vector_type& rhs, const scalar_type& s, DEG depth);
    tangent_vector& mul_scal_prod(const tangent_vector& rhs, const scalar_type& s);
    tangent_vector& mul_scal_prod(const vector_type& rhs, const scalar_type& s);

    tangent_vector& mul_scal_div(const tangent_vector& rhs, const rational_type& s, DEG depth);
    tangent_vector& mul_scal_div(const vector_type& rhs, const rational_type& s, DEG depth);
    tangent_vector& mul_scal_div(const tangent_vector& rhs, const rational_type& s);
    tangent_vector& mul_scal_div(const vector_type& rhs, const rational_type& s);
};


template <typename Coeffs, DEG Width, DEG Depth, template <typename, typename, typename...> class VectorType, typename... Args>
using tangent_tensor = tangent_vector<free_tensor<Coeffs, Width, Depth, VectorType, Args...>;



template<typename Vector, typename TangentVector>
tangent_vector<Vector, TangentVector>&
tangent_vector<Vector, TangentVector>::mul_scal_prod(const tangent_vector& rhs, const scalar_type& s, DEG depth)
{
    const auto& rhs_vec = static_cast<const vector_type&>(rhs);
    auto& this_vec = static_cast<vector_type&>(*this);
    vector_type::mul_scal_prod(rhs_vec, s, depth);
    m_tangent = this_vec.mul_scal_prod(rhs.tangent(), s, depth) + tangent().mul_scal_prod(rhs_vec, s, depth);
    return *this;
}

template<typename Vector, typename TangentVector>
tangent_vector<Vector, TangentVector>&
tangent_vector<Vector, TangentVector>::mul_scal_prod(const tangent_vector& rhs, const scalar_type& s)
{
    const auto& rhs_vec = static_cast<const vector_type&>(rhs);
    auto& this_vec = static_cast<vector_type&>(*this);
    vector_type::mul_scal_prod(rhs_vec, s);
    m_tangent = this_vec.mul_scal_prod(rhs.tangent(), s) + tangent().mul_scal_prod(rhs_vec, s);
    return *this;
}

template<typename Vector, typename TangentVector>
tangent_vector<Vector, TangentVector>&
tangent_vector<Vector, TangentVector>::mul_scal_div(const tangent_vector& rhs, const rational_type& s)
{
    const auto& rhs_vec = static_cast<const vector_type&>(rhs);
    auto& this_vec = static_cast<vector_type&>(*this);
    vector_type::mul_scal_div(rhs_vec, s);
    m_tangent = this_vec.mul_scal_div(rhs.tangent(), s) + tangent().mul_scal_div(rhs_vec, s);
    return *this;
}

template<typename Vector, typename TangentVector>
tangent_vector<Vector, TangentVector>&
tangent_vector<Vector, TangentVector>::mul_scal_div(const tangent_vector& rhs, const rational_type& s, DEG depth)
{
    const auto& rhs_vec = static_cast<const vector_type&>(rhs);
    auto& this_vec = static_cast<vector_type&>(*this);
    vector_type::mul_scal_div(rhs_vec, s, depth);
    m_tangent = this_vec.mul_scal_div(rhs.tangent(), s, depth) + tangent().mul_scal_div(rhs_vec, s, depth);
    return *this;
}

template<typename Vector, typename TangentVector>
tangent_vector<Vector, TangentVector>&
tangent_vector<Vector, TangentVector>::mul_scal_prod(const vector_type& rhs_vec, const scalar_type& s, DEG depth)
{
    auto& this_vec = static_cast<vector_type&>(*this);
    vector_type::mul_scal_prod(rhs_vec, s, depth);
    m_tangent.mul_scal_prod(rhs_vec, s, depth);
    return *this;
}
template<typename Vector, typename TangentVector>
tangent_vector<Vector, TangentVector>&
tangent_vector<Vector, TangentVector>::mul_scal_prod(const vector_type& rhs_vec, const scalar_type& s)
{
    auto& this_vec = static_cast<vector_type&>(*this);
    vector_type::mul_scal_prod(rhs_vec, s);
    m_tangent.mul_scal_prod(rhs_vec, s);
    return *this;
}
template<typename Vector, typename TangentVector>
tangent_vector<Vector, TangentVector>&
tangent_vector<Vector, TangentVector>::mul_scal_div(const vector_type& rhs_vec, const rational_type& s, DEG depth)
{
    auto& this_vec = static_cast<vector_type&>(*this);
    vector_type::mul_scal_div(rhs_vec, s, depth);
    m_tangent.mul_scal_div(rhs_vec, s, depth);
}
template<typename Vector, typename TangentVector>
tangent_vector<Vector, TangentVector>&
tangent_vector<Vector, TangentVector>::mul_scal_div(const vector_type& rhs_vec, const rational_type& s)
{
    auto& this_vec = static_cast<vector_type&>(*this);
    vector_type::mul_scal_div(rhs_vec, s);
    m_tangent.mul_scal_div(rhs_vec, s);
    return *this;
}

} // namespace alg

#endif//LIBALGEBRA_LIBALGEBRA_TANGENTS_H_
