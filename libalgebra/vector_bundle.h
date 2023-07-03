//
// Created by user on 07/07/22.
//

#ifndef LIBALGEBRA_LIBALGEBRA_VECTOR_BUNDLE_H_
#define LIBALGEBRA_LIBALGEBRA_VECTOR_BUNDLE_H_

#include "algebra.h"
#include "implementation_types.h"
#include "tensor.h"
#include "vectors.h"
#include "scalar_bundle.h"

#include <iosfwd>
#include <type_traits>
#include <utility>

#ifdef LIBALGEBRA_ENABLE_SERIALIZATION

#include <boost/serialization/serialization.hpp>

#endif

namespace alg {

namespace dtl {

template<typename Vector, typename Fibre, typename Derived>
class vector_bundle_base {
protected:
    Vector m_base;
    Fibre m_fibre;

public:
    static const typename Vector::BASIS basis;

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

    using scalar_bundle_type = scalar_bundle<coeff_type, fibre_coeff_type>;
    using rational_bundle_type = rational_bundle<coeff_type, fibre_coeff_type>;


    // Legacy declarations
    using BASIS = basis_type;
    using SCALAR = scalar_type;
    using SCA = SCALAR;
    using RATIONAL = rational_type;
    using RAT = RATIONAL;
    using KEY = key_type;
    using coefficient_field = coeff_type;

    static_assert(
            alg::utils::is_vector_type<Vector>::value,
//            std::is_base_of<vectors::vector<basis_type, coeff_type>, Vector>::value,
            "Vector must be a vector type");
    static_assert(
            alg::utils::is_vector_type<Fibre>::value,
//            std::is_base_of<vectors::vector<fibre_basis_type, fibre_coeff_type>, Fibre>::value,
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
    vector_bundle_base(const vector_bundle_base &) = default;
    vector_bundle_base(vector_bundle_base &&) noexcept = default;

    //vector_bundle_base& operator=(const vector_bundle_base& other) = default;
    vector_bundle_base& operator=(vector_bundle_base& other) noexcept = default;


    template<typename... Args, typename DummyArg = std::enable_if_t<(sizeof...(Args)>1)>>
    explicit vector_bundle_base(Args&&... args) : m_base(std::forward<Args>(args)...), m_fibre()
    {}

    template<typename Arg, typename DummyArg = std::enable_if_t<!(std::is_same<Arg, vector_bundle_base>::value)>>
    explicit vector_bundle_base(Arg&& arg) : m_base(std::forward<Arg>(arg)), m_fibre()
    {}

    explicit vector_bundle_base(const Vector &point) : m_base(point), m_fibre() {}

    explicit vector_bundle_base(Vector &&point) : m_base(std::move(point)), m_fibre() {}

    vector_bundle_base(Vector &&point, Fibre &&tangent)
            : m_base(std::move(point)), m_fibre(std::move(tangent)) {}

    vector_bundle_base(Vector &&point, const Fibre &fibre)
            : m_base(std::move(point)), m_fibre(fibre) {}

    vector_bundle_base(const Vector &point, const Fibre &tangent)
            : m_base(point), m_fibre(tangent) {}

    template<typename Key, typename = typename std::enable_if<
            std::is_same<Key, key_type>::value && std::is_same<Key, fibre_key_type>::value>>
    explicit vector_bundle_base(const Key& k)
            : m_base(k), m_fibre(k) {}

    template<typename... Args>
    explicit vector_bundle_base(Vector &&point, Args &&... args)
            : m_base(std::move(point)), m_fibre(std::forward<Args>(args)...) {}

    template<typename... Args>
    explicit vector_bundle_base(const Vector &point, Args &&... args)
            : m_base(point), m_fibre(std::forward<Args>(args)...) {}

    Vector &base() noexcept { return m_base; }

    const Vector &base() const noexcept { return m_base; }

    Fibre &fibre() noexcept {
        return m_fibre;
    }

    const Fibre &fibre() const noexcept {
        return m_fibre;
    }

private:
#ifdef LIBALGEBRA_ENABLE_SERIALIZATION

    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive &ar, unsigned int const /*version*/) {
        ar & m_base;
        ar & m_fibre;
    }

#endif
public:

    void swap(vector_bundle_base& other) noexcept {
        std::swap(m_base, other.m_base);
        std::swap(m_fibre, other.m_fibre);
    }

    const scalar_type& operator[](const key_type& key) const
    { return m_base[key]; }

    Derived& add_scal_prod(const typename Vector::KEY& key, const scalar_type& s);
    Derived& sub_scal_prod(const typename Vector::KEY& key, const scalar_type& s);

    Derived &add_scal_prod(const vector_bundle_base &rhs, const scalar_type &s);
    Derived &sub_scal_prod(const vector_bundle_base &rhs, const scalar_type &s);
    Derived &add_scal_div(const vector_bundle_base &rhs, const rational_type &s);
    Derived &sub_scal_div(const vector_bundle_base &rhs, const rational_type &s);
    Derived &add_scal_prod(const Vector &rhs, const scalar_type &s);
    Derived &sub_scal_prod(const Vector &rhs, const scalar_type &s);
    Derived &add_scal_div(const Vector &rhs, const rational_type &s);
    Derived &sub_scal_div(const Vector &rhs, const rational_type &s);

    template<typename OtherVector>
    typename std::enable_if<!std::is_base_of<vector_bundle_base<Vector, Fibre, Derived>, OtherVector>::value, Derived>::type &
    add_scal_prod(const OtherVector &rhs, const scalar_type &s);

    Derived &mul_scal_prod(const vector_bundle_base &rhs, const scalar_type &s, DEG depth);

    Derived &mul_scal_prod(const vector_type &rhs, const scalar_type &s, DEG depth);

    Derived &mul_scal_prod(const vector_bundle_base &rhs, const scalar_type &s);

    Derived &mul_scal_prod(const vector_type &rhs, const scalar_type &s);

    Derived &mul_scal_div(const vector_bundle_base &rhs, const rational_type &s, DEG depth);

    Derived &mul_scal_div(const vector_type &rhs, const rational_type &s, DEG depth);

    Derived &mul_scal_div(const vector_bundle_base &rhs, const rational_type &s);

    Derived &mul_scal_div(const vector_type &rhs, const rational_type &s);

    Derived& add_mul(const vector_bundle_base& lhs, const vector_bundle_base& rhs);
    Derived& sub_mul(const vector_bundle_base& lhs, const vector_bundle_base& rhs);


    scalar_type NormL1() const { return m_base.NormL1(); }

};


template <typename Vector, typename Fibre, typename Derived>
const typename Vector::BASIS vector_bundle_base<Vector, Fibre, Derived>::basis = Vector::basis;

template<typename Vector, typename Fibre, typename Derived>
Derived operator-(const vector_bundle_base<Vector, Fibre, Derived>& arg)
{
    return {-arg.base(), -arg.fibre()};
}

template<typename LBase, typename LFibre, typename LDerived,
         typename RBase, typename RFibre, typename RDerived>
LDerived operator+(const vector_bundle_base<LBase, LFibre, LDerived>& lhs,
                   const vector_bundle_base<RBase, RFibre, RDerived>& rhs)
{
    return {lhs.base() + rhs.base(), lhs.fibre() + rhs.fibre()};
}

template<typename Vector, typename Fibre, typename Derived>
Derived operator+(const vector_bundle_base<Vector, Fibre, Derived>& lhs, const Vector& rhs)
{
    return {lhs.base() + rhs, lhs.fibre()};
}

template<typename Vector, typename Fibre, typename Derived>
Derived operator+(const Vector& lhs, const vector_bundle_base<Vector, Fibre, Derived>& rhs)
{
    return {lhs + rhs.base(), rhs.fibre()};
}

template<typename LBase, typename LFibre, typename LDerived,
         typename RBase, typename RFibre, typename RDerived>
LDerived operator-(const vector_bundle_base<LBase, LFibre, LDerived>& lhs,
                   const vector_bundle_base<RBase, RFibre, RDerived>& rhs)
{
    return {lhs.base() - rhs.base(), lhs.fibre() - rhs.fibre()};
}

template<typename Vector, typename Fibre, typename Derived>
Derived operator-(const vector_bundle_base<Vector, Fibre, Derived>& lhs, const Vector& rhs)
{
    return {lhs.base() - rhs, lhs.fibre()};
}

template<typename Vector, typename Fibre, typename Derived>
Derived operator-(const Vector& lhs, const vector_bundle_base<Vector, Fibre, Derived>& rhs)
{
    return {lhs - rhs.base(), -rhs.fibre()};
}

/*
 * For scalar multiplication, the vector and fibre types might have different scalar types, so we define
 * scalar multiplication with respect to both and let the compiler figure out which one is more permissive.
 * The same is true for rational division.
 */

template<typename Vector, typename Fibre, typename Derived>
Derived
operator*(const vector_bundle_base<Vector, Fibre, Derived>& lhs, const typename Vector::SCALAR& rhs)
{
    return {lhs.base() * rhs, lhs.fibre() * rhs};
}

template<typename Vector, typename Fibre, typename Derived>
Derived operator*(const typename Vector::SCALAR& lhs, const vector_bundle_base<Vector, Fibre, Derived>& rhs)
{
    return {lhs * rhs.base(), lhs * rhs.fibre()};
}

template<typename Vector, typename Fibre, typename Derived>
Derived operator/(const vector_bundle_base<Vector, Fibre, Derived>& lhs, const typename Vector::RATIONAL& rhs)
{
    return {lhs.base() / rhs, lhs.fibre() / rhs};
}

template<typename LBase, typename LFibre, typename LDerived,
         typename RBase, typename RFibre, typename RDerived>
LDerived operator*(
        const vector_bundle_base<LBase, LFibre, LDerived>& left,
        const vector_bundle_base<RBase, RFibre, RDerived>& right)
{
    const auto& lhs_vec = left.base();
    const auto& rhs_vec = right.base();
    return {lhs_vec * rhs_vec, lhs_vec * right.fibre() + left.fibre() * rhs_vec};
}

template<typename Vector, typename Fibre, typename Derived>
Derived operator*(
        const vector_bundle_base<Vector, Fibre, Derived>& left,
        const Vector& rhs_vec)
{
    const auto& lhs_vec = left.base();
    return {lhs_vec * rhs_vec, left.fibre() * rhs_vec};
}

template<typename Vector, typename Fibre, typename Derived>
Derived operator*(
        const Vector& lhs_vec,
        const vector_bundle_base<Vector, Fibre, Derived>& right)
{
    const auto& rhs_vec = right.base();
    return {lhs_vec * rhs_vec, lhs_vec * right.fibre()};
}

/*
 * Now all the in-place operations.
 */

template<typename LBase, typename LFibre, typename LDerived,
         typename RBase, typename RFibre, typename RDerived>
LDerived&
operator+=(vector_bundle_base<LBase, LFibre, LDerived>& lhs,
           const vector_bundle_base<RBase, RFibre, RDerived>& rhs)
{
    lhs.base() += rhs.base();
    lhs.fibre() += rhs.fibre();
    return static_cast<LDerived&>(lhs);
}

template<typename Vector, typename Fibre, typename Derived>
Derived& operator+=(vector_bundle_base<Vector, Fibre, Derived>& lhs, const Vector& rhs)
{
    lhs.base() += rhs;
    return static_cast<Derived&>(lhs);
}

template<typename LBase, typename LFibre, typename LDerived,
         typename RBase, typename RFibre, typename RDerived>
LDerived&
operator-=(vector_bundle_base<LBase, LFibre, LDerived>& lhs,
           const vector_bundle_base<RBase, RFibre, RDerived>& rhs)
{
    lhs.base() -= rhs.base();
    lhs.fibre() -= rhs.fibre();
    return static_cast<LDerived&>(lhs);
}

template<typename Vector, typename Fibre, typename Derived>
Derived& operator-=(vector_bundle_base<Vector, Fibre, Derived>& lhs, const Vector& rhs)
{
    lhs.base() -= rhs;
    return static_cast<Derived&>(lhs);
}

template<typename Vector, typename Fibre, typename Derived>
Derived& operator*=(vector_bundle_base<Vector, Fibre, Derived>& lhs, const typename Vector::SCALAR& rhs)
{
    lhs.base() *= rhs;
    lhs.fibre() *= rhs;
    return static_cast<Derived&>(lhs);
}

template<typename Vector, typename Fibre, typename Derived>
Derived& operator/=(vector_bundle_base<Vector, Fibre, Derived>& lhs, const typename Vector::RATIONAL& rhs)
{
    lhs.base() /= rhs;
    lhs.fibre() /= rhs;
    return static_cast<Derived&>(lhs);
}

template<typename LBase, typename LFibre, typename LDerived,
         typename RBase, typename RFibre, typename RDerived>
LDerived& operator*=(
        vector_bundle_base<LBase, LFibre, LDerived>& lhs,
        const vector_bundle_base<RBase, RFibre, RDerived>& rhs)
{
    auto& lhs_vec = lhs.base();
    const auto& rhs_vec = rhs.base();

    lhs.fibre() *= rhs_vec;
    lhs.fibre() += lhs_vec * rhs.fibre();

    // Do this operation last because otherwise it messes with the
    // calculation of the fibre.
    lhs_vec *= rhs_vec;
    return static_cast<LDerived&>(lhs);
}

template<typename Vector, typename Fibre, typename Derived>
Derived& operator*=(
        vector_bundle_base<Vector, Fibre, Derived>& lhs,
        const Vector& rhs_vec)
{
    auto& lhs_vec = lhs.base();

    lhs_vec *= rhs_vec;
    lhs.fibre() *= rhs_vec;
    return static_cast<Derived&>(lhs);
}

template<typename LBase, typename LFibre, typename LDerived,
         typename RBase, typename RFibre, typename RDerived>
bool operator==(const vector_bundle_base<LBase, LFibre, LDerived>& lhs,
                const vector_bundle_base<RBase, RFibre, RDerived>& rhs)
{
    return (lhs.base() == rhs.base()) && (lhs.fibre() == rhs.fibre());
}

template<typename LBase, typename LFibre, typename LDerived,
         typename RBase, typename RFibre, typename RDerived>
bool operator!=(const vector_bundle_base<LBase, LFibre, LDerived>& lhs,
                const vector_bundle_base<RBase, RFibre, RDerived>& rhs)
{
    return !operator==(lhs, rhs);
}

template<typename Vector, typename Fibre, typename Derived>
std::ostream& operator<<(std::ostream& os, const vector_bundle_base<Vector, Fibre, Derived>& arg)
{
    return os << '('
              << arg.base()
              << ", "
              << arg.fibre()
              << ')';
}

template<typename LBase, typename LFibre, typename LDerived,
         typename RBase, typename RFibre, typename RDerived>
typename std::enable_if<is_algebra<LBase>(), LDerived>::type
commutator(const vector_bundle_base<LBase, LFibre, LDerived>& x,
           const vector_bundle_base<RBase, RFibre, RDerived>& y)
{
    LDerived result(x * y);
    result.add_mul(y, x);
    return result;
}

template <typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::add_scal_prod(const typename Vector::KEY& key, const scalar_type& s) {
    m_base.add_scal_prod(key, s);
    return static_cast<Derived&>(*this);
}
template <typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::sub_scal_prod(const typename Vector::KEY& key, const scalar_type& s) {
    m_base.sub_scal_prod(key, s);
    return static_cast<Derived&>(*this);
}

template<typename Vector, typename Fibre, typename Derived>
Derived&
vector_bundle_base<Vector, Fibre, Derived>::add_scal_prod(const vector_bundle_base& rhs, const scalar_type& s)
{
    m_base.add_scal_prod(rhs.base(), s);
    m_fibre.add_scal_prod(rhs.m_fibre, s);
    return static_cast<Derived&>(*this);
}

template<typename Vector, typename Fibre, typename Derived>
Derived&
vector_bundle_base<Vector, Fibre, Derived>::sub_scal_prod(const vector_bundle_base& rhs, const scalar_type& s)
{
    m_base.sub_scal_prod(rhs.base(), s);
    m_fibre.sub_scal_prod(rhs.m_fibre, s);
    return static_cast<Derived&>(*this);
}

template<typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::add_scal_div(const vector_bundle_base& rhs,
                                                                  const rational_type& s)
{
    m_base.add_scal_div(rhs.base(), s);
    m_fibre.add_scal_div(rhs.m_fibre, s);
    return static_cast<Derived&>(*this);
}

template<typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::sub_scal_div(const vector_bundle_base& rhs,
                                                                  const rational_type& s)
{
    m_base.sub_scal_div(rhs.base(), s);
    m_fibre.sub_scal_div(rhs.m_fibre, s);
    return static_cast<Derived&>(*this);
}

template<typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::add_scal_prod(const Vector& rhs, const scalar_type& s)
{
    m_base.add_scal_prod(rhs, s);
    return static_cast<Derived&>(*this);
}

template<typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::sub_scal_prod(const Vector& rhs, const scalar_type& s)
{
    m_base.sub_scal_prod(rhs, s);
    return static_cast<Derived&>(*this);
}

template<typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::add_scal_div(const Vector& rhs, const rational_type& s)
{
    m_base.add_scal_div(rhs, s);
    return static_cast<Derived&>(*this);
}

template<typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::sub_scal_div(const Vector& rhs, const rational_type& s)
{
    m_base.sub_scal_div(rhs, s);
    return static_cast<Derived&>(*this);
}

template<typename Vector, typename Fibre, typename Derived>
template<typename OtherVector>
typename std::enable_if<!std::is_base_of<vector_bundle_base<Vector, Fibre, Derived>, OtherVector>::value, Derived>::type&
vector_bundle_base<Vector, Fibre, Derived>::add_scal_prod(const OtherVector& rhs, const scalar_type& s)
{
    m_base.add_scal_prod(rhs, s);
    return static_cast<Derived&>(*this);
}

template<typename Vector, typename Fibre, typename Derived>
Derived&
vector_bundle_base<Vector, Fibre, Derived>::mul_scal_prod(const vector_bundle_base& rhs, const scalar_type& s,
                                                          DEG depth)
{
    m_fibre.mul_scal_prod(rhs.base(), s, depth);
    m_fibre += m_base * (rhs.m_fibre * s);
    m_base.mul_scal_prod(rhs.base(), s, depth);
    return static_cast<Derived&>(*this);
}

template<typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::mul_scal_prod(const vector_type& rhs, const scalar_type& s,
                                                                   DEG depth)
{
    m_fibre.mul_scal_prod(rhs, s, depth);
    m_base.mul_scal_prod(rhs, s, depth);
    return static_cast<Derived&>(*this);
}

template<typename Vector, typename Fibre, typename Derived>
Derived&
vector_bundle_base<Vector, Fibre, Derived>::mul_scal_prod(const vector_bundle_base& rhs, const scalar_type& s)
{
    m_fibre.mul_scal_prod(rhs.base(), s);
    m_fibre += m_base * (rhs.m_fibre * s);
    m_base.mul_scal_prod(rhs.base(), s);
    return static_cast<Derived&>(*this);
}

template<typename Vector, typename Fibre, typename Derived>
Derived&
vector_bundle_base<Vector, Fibre, Derived>::mul_scal_prod(const vector_type& rhs, const scalar_type& s)
{
    m_base.mul_scal_prod(rhs, s);
    m_fibre.mul_scal_prod(rhs, s);
    return static_cast<Derived&>(*this);
}

template<typename Vector, typename Fibre, typename Derived>
Derived&
vector_bundle_base<Vector, Fibre, Derived>::mul_scal_div(const vector_bundle_base& rhs, const rational_type& s,
                                                         DEG depth)
{
    m_fibre.mul_scal_div(rhs.base(), s, depth);
    m_fibre += m_base * (rhs.m_fibre / s);
    m_base.mul_scal_div(rhs.base(), s, depth);
    return static_cast<Derived&>(*this);
}

template<typename Vector, typename Fibre, typename Derived>
Derived&
vector_bundle_base<Vector, Fibre, Derived>::mul_scal_div(const vector_type& rhs, const rational_type& s,
                                                         DEG depth)
{
    m_base.mul_scal_div(rhs, s, depth);
    m_fibre.mul_scal_div(rhs, s, depth);
    return static_cast<Derived&>(*this);
}

template<typename Vector, typename Fibre, typename Derived>
Derived& vector_bundle_base<Vector, Fibre, Derived>::mul_scal_div(const vector_bundle_base& rhs,
                                                                  const rational_type& s)
{
    m_fibre.mul_scal_div(rhs.base(), s);
    m_fibre += m_base * (rhs.m_fibre / s);
    m_base.mul_scal_div(rhs.base(), s);
    return static_cast<Derived&>(*this);
}

template<typename Vector, typename Fibre, typename Derived>
Derived&
vector_bundle_base<Vector, Fibre, Derived>::mul_scal_div(const vector_type& rhs, const rational_type& s)
{
    m_base.mul_scal_div(rhs, s);
    m_fibre.mul_scal_div(rhs, s);
    return static_cast<Derived&>(*this);
}

template<typename Vector, typename Fibre, typename Derived>
Derived&
vector_bundle_base<Vector, Fibre, Derived>::add_mul(const vector_bundle_base& lhs, const vector_bundle_base& rhs)
{
    // Operation is z = z + x*y
    // derivative is z' = z' + (x*y' + x'*y)
    m_fibre += lhs.base() * rhs.fibre() + lhs.fibre() * rhs.base();
    m_base.sub_mul(lhs.base(), rhs.base());

    return static_cast<Derived&>(*this);
}

template<typename Vector, typename Fibre, typename Derived>
Derived&
vector_bundle_base<Vector, Fibre, Derived>::sub_mul(const vector_bundle_base& lhs, const vector_bundle_base& rhs)
{
    // Operation is z = z - x*y
    // derivative is z' = z' - (x*y' + x'*y)
    m_fibre -= lhs.base() * rhs.fibre() + lhs.fibre() * rhs.base();
    m_base.sub_mul(lhs.base(), rhs.base());

    return static_cast<Derived&>(*this);
}

}// namespace dtl






template<typename Vector, typename Fibre = Vector>
class vector_bundle : public dtl::vector_bundle_base<Vector, Fibre, vector_bundle<Vector, Fibre>> {
public:
    using bundle_base = dtl::vector_bundle_base<Vector, Fibre, vector_bundle>;

    using bundle_base::bundle_base;

    using typename bundle_base::vector_type;
    using typename bundle_base::fibre_vector_type;

    using typename bundle_base::basis_type;
    using typename bundle_base::key_type;
    using typename bundle_base::coeff_type;
    using typename bundle_base::scalar_type;
    using typename bundle_base::rational_type;

    using typename bundle_base::fibre_basis_type;
    using typename bundle_base::fibre_key_type;
    using typename bundle_base::fibre_coeff_type;
    using typename bundle_base::fibre_scalar_type;
    using typename bundle_base::fibre_rational_type;

    using typename bundle_base::BASIS;
    using typename bundle_base::KEY;
    using typename bundle_base::SCALAR;
    using typename bundle_base::RATIONAL;
    using typename bundle_base::coefficient_field;
private:
#ifdef LIBALGEBRA_ENABLE_SERIALIZATION

    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive &ar, unsigned int const /*version*/) {
        ar & boost::serialization::base_object<bundle_base>(*this);
    }

#endif
};




template<typename Coeffs, DEG Width, DEG Depth, template<typename, typename, typename...> class VectorType, template <DEG, DEG> class FTMultiplication, typename... Args>
class vector_bundle<free_tensor<Coeffs, Width, Depth, VectorType, FTMultiplication, Args...>,
                    free_tensor<Coeffs, Width, Depth, VectorType, FTMultiplication, Args...>>
        : public dtl::vector_bundle_base<free_tensor<Coeffs, Width, Depth,  VectorType, FTMultiplication, Args...>,
                                         free_tensor<Coeffs, Width, Depth, VectorType, FTMultiplication, Args...>,
                                         vector_bundle<
                                                free_tensor < Coeffs, Width, Depth, VectorType, FTMultiplication, Args...>,
                                                free_tensor<Coeffs, Width, Depth, VectorType, FTMultiplication, Args...>
                                            >
                                         >
{
#ifdef LIBALGEBRA_ENABLE_SERIALIZATION

    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive &ar, unsigned int const /*version*/) {
        ar & boost::serialization::base_object<bundle_base>(*this);
    }

#endif


    using vector_t = free_tensor<Coeffs, Width, Depth, VectorType, FTMultiplication, Args...>;
    using fibre_t = vector_t;
public:

    using bundle_base = dtl::vector_bundle_base<vector_t, fibre_t, vector_bundle>;

    using bundle_base::bundle_base;

    using typename bundle_base::vector_type;
    using typename bundle_base::fibre_vector_type;

    using typename bundle_base::basis_type;
    using typename bundle_base::key_type;
    using typename bundle_base::coeff_type;
    using typename bundle_base::scalar_type;
    using typename bundle_base::rational_type;

    using typename bundle_base::fibre_basis_type;
    using typename bundle_base::fibre_key_type;
    using typename bundle_base::fibre_coeff_type;
    using typename bundle_base::fibre_scalar_type;
    using typename bundle_base::fibre_rational_type;

    using typename bundle_base::BASIS;
    using typename bundle_base::KEY;
    using typename bundle_base::SCALAR;
    using typename bundle_base::RATIONAL;
    using typename bundle_base::coefficient_field;

private:


    template<typename V>
    static void resize_for_degree_impl(V &, DEG) {}

    static void resize_for_degree_impl(vectors::dense_vector<basis_type, coeff_type> &vec, DEG degree) {
        vec.resize_to_degree(degree);
    }

    static void resize_for_degree(vector_bundle &bundle, DEG depth) {
        resize_for_degree_impl(bundle.base().base_vector(), depth);
        resize_for_degree_impl(bundle.fibre().base_vector(), depth);
    }

public:
    vector_bundle fmexp(const vector_bundle& arg) const
    {
        vector_bundle result(*this);
        vector_bundle x(arg);
        key_type kunit;

        //        const auto& self = static_cast<const vector_type&>(*this);

        auto unit_elt = x.base().find(kunit);
        if (unit_elt != x.base().end() && unit_elt->value() != coeff_type::zero) {
            x.base().erase(unit_elt);
        }

        for (DEG i = Depth; i >= 1; --i) {
            result.mul_scal_div(x, rational_type(i), Depth - i + 1);
            result += *this;
        }

        return result;
    }

    vector_bundle fmexp(const vector_type& arg) const
    {
        vector_bundle result(*this);
        vector_type x(arg);

        key_type kunit;
        auto unit_elt = x.find(kunit);
        if (unit_elt != x.end() && unit_elt->value() != coeff_type::zero) {
            x.erase(unit_elt);
        }

        for (DEG i = Depth; i >= 1; --i) {
            result.mul_scal_div(x, rational_type(i), Depth - i + 1);
            result += *this;
        }

        return result;
    }

    vector_bundle& fmexp_inplace(const vector_bundle& arg)
    {
        vector_bundle self(*this);
        vector_bundle x(arg);

        key_type kunit;

        auto unit_elt = x.base().find(kunit);
        if (unit_elt != x.base().end() && unit_elt->value() != coeff_type::zero) {
            x.base().erase(unit_elt);
        }

        for (DEG i = Depth; i >= 1; --i) {
            bundle_base::mul_scal_div(x, rational_type(i), Depth - i + 1);
            *this += self;
        }

        return *this;
    }

    vector_bundle& fmexp_inplace(const vector_type& arg)
    {
        vector_bundle self(*this);
        vector_type x(arg);

        key_type kunit;
        auto unit_elt = x.find(kunit);
        if (unit_elt != x.end() && unit_elt->value() != coeff_type::zero) {
            x.erase(unit_elt);
        }

        for (DEG i = Depth; i >= 1; --i) {
            bundle_base::mul_scal_div(x, rational_type(i), Depth - i + 1);
            *this += self;
        }

        return *this;
    }

    friend vector_bundle antipode(const vector_bundle& arg)
    {
        return {antipode(arg.base()), antipode(arg.fibre())};
    }

    friend vector_bundle exp(const vector_bundle& arg) {
        // Computes the truncated exponential of arg
        // 1 + arg + arg^2/2! + ... + arg^n/n! where n = max_degree
        KEY kunit;
        vector_bundle tunit(kunit);
        vector_bundle result(tunit);

        resize_for_degree(result, Depth);

        for (DEG i = Depth; i >= 1; --i) {
            result.mul_scal_div(arg, rational_type(i));
            result += tunit;
        }
        return result;
    }


    friend vector_bundle log(const vector_bundle& arg)
    {
        // Computes the truncated log of arg up to degree max_degree
        // The coef. of the constant term (empty word in the monoid) of arg
        // is forced to 1.
        // log(arg) = log(1+x) = x - x^2/2 + ... + (-1)^(n+1) x^n/n.
        // max_degree must be > 0

        KEY kunit;
        vector_bundle tunit(kunit);
        vector_bundle x(arg);
        auto it = x.base().find(kunit);
        if (it != x.base().end()) {
            x.base().erase(it);
        }
        vector_bundle result;

        for (DEG i = Depth; i >= 1; --i) {
            if (i % 2 == 0) {
                result.sub_scal_div(tunit, typename Coeffs::Q(i));
            }
            else {
                result.add_scal_div(tunit, typename Coeffs::Q(i));
            }
            result *= x;
        }

        return result;
    }


    friend vector_bundle inverse(const vector_bundle &arg) {
        // Computes the truncated inverse of arg up to degree max_degree
        // An exception is thrown if the leading term is zero.
        // the module assumes
        // (a+x)^(-1) = (a(1+x/a))^(-1)
        //  = a^(-1)(1 - x/a + x^2/a^2 + ... + (-1)^(n) x^n/a^n)
        // = a^(-1) - x/a*[a^(-1)(1 - x/a + x^2/a^2 + ... + (-1)^(n)
        // x^(n-1)/a^(n-1)))]. S_n = a^(-1) + z S_{n-1}; z = - x/a ; S_0 = a^(-1)
        // max_degree must be > 0

        KEY kunit;
        scalar_type a(0);
        vector_bundle x;
        vector_bundle z(a);

        auto it = arg.base().find(kunit);
        if (it == arg.base().end()) {
            // const term a is 0;
            throw std::invalid_argument("divide-by-zero");
        } else {
            a = (*it).value();
            x = arg;
            x.base().erase(kunit);
        }

        // S_n = a + z S_{ n - 1 }; z = -x / a; S_0 = a
        //
        //  the nonzero scalar component a of the tensor arg restored to a tensor
        vector_bundle free_tensor_a_inverse(scalar_type(1) / a);
        vector_bundle result(free_tensor_a_inverse);
        // z := - x/a
        z.sub_scal_div(x, a);
        // the iteration
        for (DEG i = 0; i != Depth; ++i) {
            auto tmp = z * result;
            result = free_tensor_a_inverse + z * result;
        }
        return result;
    }

};



namespace dtl {

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
