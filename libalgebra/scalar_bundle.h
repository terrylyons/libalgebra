//
// Created by user on 05/05/23.
//

#ifndef LIBALGEBRA_LIBALGEBRA_SCALAR_BUNDLE_H_
#define LIBALGEBRA_LIBALGEBRA_SCALAR_BUNDLE_H_

#include "implementation_types.h"

#include <iosfwd>
#include <type_traits>
#include <utility>

#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/serialization.hpp>
#endif

namespace alg {
namespace dtl {

template<typename Base, typename Fibre, typename Derived>
class scalar_bundle_base
{
    Base m_base;
    Fibre m_fibre;

public:
    using base_type = Base;
    using FibreType = Fibre;

    constexpr scalar_bundle_base(Base&& base, Fibre&& fibre)
        : m_base(std::move(base)), m_fibre(std::move(fibre))
    {}

    explicit constexpr scalar_bundle_base(Base&& base) : m_base(std::move(base)), m_fibre()
    {}

    template <typename... Args>
    constexpr explicit scalar_bundle_base(Args&&... args)
            : m_base(std::forward<Args>(args)...), m_fibre()
    {}

    template <typename Arg>
    constexpr scalar_bundle_base& operator=(const Arg& other) {
        m_base = static_cast<Base>(other);
        return *this;
    }

    constexpr Base& base() noexcept { return m_base; }
    constexpr Fibre& fibre() noexcept { return m_fibre; }
    constexpr const Base& base() const noexcept { return m_base; }
    constexpr const Fibre& fibre() const noexcept { return m_fibre; }

    void swap(scalar_bundle_base& other)
    {
        std::swap(m_base, other.m_base);
        std::swap(m_fibre, other.m_fibre);
    }

private:
#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive& ar, unsigned int const /*version*/)
    {
        ar & m_base;
        ar & m_fibre;
    }
#endif
};

}// namespace dtl

template<typename BaseRing, typename FibreRing>
class scalar_bundle
    : public dtl::scalar_bundle_base<typename BaseRing::S,
                                     typename FibreRing::S,
                                     scalar_bundle<BaseRing, FibreRing>>
{
    using bundle_base = dtl::scalar_bundle_base<typename BaseRing::S,
                                                typename FibreRing::S,
                                                scalar_bundle>;

public:

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

template<typename BaseRing, typename FibreRing>
class rational_bundle
    : public dtl::scalar_bundle_base<typename BaseRing::Q,
                                     typename FibreRing::Q,
                                     scalar_bundle<BaseRing, FibreRing>>
{
    using bundle_base = dtl::scalar_bundle_base<typename BaseRing::Q,
                                                typename FibreRing::Q,
                                                rational_bundle>;

public:

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





// Implementations follow
namespace dtl {

template<typename B, typename F, typename Derived>
constexpr Derived operator-(const scalar_bundle_base<B, F, Derived>& arg)
{
    return {-arg.base(), -arg.fibre()};
}

template<typename B, typename F, typename Derived>
constexpr Derived operator+(const scalar_bundle_base<B, F, Derived>& lhs,
                            const scalar_bundle_base<B, F, Derived>& rhs) {
    return { lhs.base() + rhs.base(), lhs.fibre() + rhs.fibre() };
}

template<typename B, typename F, typename Derived>
constexpr Derived operator+(const scalar_bundle_base<B, F, Derived>& lhs,
                            const B& rhs) {
    return { lhs.base() + rhs, lhs.fibre() };
}

template<typename B, typename F, typename Derived>
constexpr Derived operator+(const B& lhs,
                            const scalar_bundle_base<B, F, Derived>& rhs) {
    return { lhs + rhs.base(), rhs.fibre() };
}

template<typename B, typename F, typename Derived>
constexpr Derived operator-(const scalar_bundle_base<B, F, Derived>& lhs,
                            const scalar_bundle_base<B, F, Derived>& rhs) {
    return { lhs.base() - rhs.base(), lhs.fibre() - rhs.fibre() };
}

template<typename B, typename F, typename Derived>
constexpr Derived operator-(const scalar_bundle_base<B, F, Derived>& lhs,
                            const B& rhs)
{
    return {lhs.base() - rhs, lhs.fibre()};
}

template<typename B, typename F, typename Derived>
constexpr Derived operator-(const B& lhs,
                            const scalar_bundle_base<B, F, Derived>& rhs)
{
    return {lhs - rhs.base(), rhs.fibre()};
}

template<typename B, typename F, typename Derived>
constexpr Derived operator*(const scalar_bundle_base<B, F, Derived>& lhs,
                            const scalar_bundle_base<B, F, Derived>& rhs) {
    return { lhs.base() * rhs.base(), lhs.fibre() * rhs.base() + lhs.base() * rhs.fibre() };
}

template<typename B, typename F, typename Derived>
constexpr Derived operator*(const scalar_bundle_base<B, F, Derived>& lhs,
                            const B& rhs)
{
    return {lhs.base() * rhs, lhs.fibre()*rhs };
}

template<typename B, typename F, typename Derived>
constexpr Derived operator*(const B& lhs,
                            const scalar_bundle_base<B, F, Derived>& rhs)
{
    return {lhs * rhs.base(), lhs*rhs.fibre()};
}


template<typename B, typename F, typename Derived>
constexpr Derived& operator+=(scalar_bundle_base<B, F, Derived>& lhs,
                              const scalar_bundle_base<B, F, Derived>& rhs) {
    lhs.base() += rhs.base();
    lhs.fibre() += rhs.fibre();
    return static_cast<Derived&>(lhs);
}

template<typename B, typename F, typename Derived>
constexpr Derived& operator+=(scalar_bundle_base<B, F, Derived>& lhs,
                              const B& rhs) {
    lhs.base() += rhs;
    return static_cast<Derived&>(lhs);
}

template<typename B, typename F, typename Derived>
constexpr Derived& operator-=(scalar_bundle_base<B, F, Derived>& lhs,
                              const scalar_bundle_base<B, F, Derived>& rhs) {
    lhs.base() -= rhs.base();
    lhs.fibre() -= rhs.fibre();
    return static_cast<Derived&>(lhs);
}

template<typename B, typename F, typename Derived>
constexpr Derived& operator-=(scalar_bundle_base<B, F, Derived>& lhs,
                              const B& rhs) {
    lhs.base() -= rhs;
    return static_cast<Derived&>(lhs);
}

template<typename B, typename F, typename Derived>
constexpr Derived& operator*=(scalar_bundle_base<B, F, Derived>& lhs,
                              const scalar_bundle_base<B, F, Derived>& rhs) {
    lhs.fibre() = lhs.base()*rhs.fibre() + lhs.fibre()*rhs.base();
    lhs.base() *= rhs.base();
    return static_cast<Derived&>(lhs);
}

template<typename B, typename F, typename Derived>
constexpr Derived& operator*=(scalar_bundle_base<B, F, Derived>& lhs,
                              const B& rhs) {
    lhs.base() *= rhs;
    lhs.fibre() *= rhs;
    return static_cast<Derived&>(lhs);
}

}// namespace dtl

template <typename BR, typename FR>
constexpr scalar_bundle<BR, FR> operator/(const scalar_bundle<BR, FR>& num,
                                          const rational_bundle<BR, FR>& den) {
    return { num.base() / den.base(),
            (num.fibre()*den.base() - num.base()*den.fibre()) / (den.fibre() * den.fibre())};
}

template <typename BR, typename FR>
constexpr scalar_bundle<BR, FR>& operator/=(scalar_bundle<BR, FR>& num,
                                            const rational_bundle<BR, FR>& den) {
    num.fibre() = (num.fibre()*den.base() - num.base()*den.fibre()) / (den.fibre()*den.fibre());
    num.base() /= den.base();
    return num;
}


}// namespace alg

#endif//LIBALGEBRA_LIBALGEBRA_SCALAR_BUNDLE_H_
