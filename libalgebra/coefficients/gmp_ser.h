//
// Created by sam on 26/10/2021.
//

#ifndef GMP_SER_GMP_SER_H
#define GMP_SER_GMP_SER_H

#include <gmp.h>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/array_wrapper.hpp>

#include <boost/multiprecision/gmp.hpp>

BOOST_SERIALIZATION_SPLIT_FREE(boost::multiprecision::backends::gmp_int)
BOOST_SERIALIZATION_SPLIT_FREE(boost::multiprecision::backends::gmp_rational)

/*
 * So the helper functions mpz_limbs_read, mpz_limbs_write, and mpz_limbs_finish were added
 * in GMP v6.0 so we have to do something different for GMP < 6.0. Reading is fine because
 * we can just grab the pointer from the struct, but writing is more tricky because it
 * involves reallocating the mpz/mpq. To avoid repeating code, we define some helpers
 * here.
 */

namespace gmp_ser_helpers {

template<typename Archive>
inline void save_to_archive(Archive& ar, mpz_srcptr mpz)
{
    mp_size_t sz = mpz_size(mpz);
    int sign = mpz_sgn(mpz);
    ar << boost::serialization::make_nvp("size", sz);
    ar << boost::serialization::make_nvp("sign", sign);

    // We're going to store this as an array of limb_t (= unsigned long int on many systems).
    // Obviously this isn't necessarily portable.
    // These are essentially equivalent, but the first is making use of mpz_limbs_read
    // for forward safety
#if __GNU_MP_VERSION >= 6
    ar << boost::serialization::make_array(mpz_limbs_read(mpz), sz);
#else
    ar << boost::serialization::make_array(mpz->_mp_d)
#endif
}

template<typename Archive>
inline void load_from_archive(Archive& ar, mpz_ptr mpz)
{
    // Expect an initialised mpz.
    mp_size_t sz = 0;
    int sign = 0;
    ar >> boost::serialization::make_nvp("size", sz);
    ar >> boost::serialization::make_nvp("sign", sign);

    if (sz > 0) {
        // We need to get write access to the limbs of mpz, which possibly means
        // allocating new space for the limbs. The better option is to use the
        // utility functions mpz_limbs_write and mpz_limbs_finish for version 6.
        // For versions <6, we have to do this by hand.
#if __GNU_MP_VERSION >= 6
        // The write function reallocates if necessary
        mp_limb_t* data = mpz_limbs_write(mpz, sz);
        // Read directly into the limb array
        ar >> boost::serialization::make_array(data, sz);
        mpz_limbs_finish(mpz, static_cast<mp_size_t>(sign) * sz);
#else
        mp_size_t signed_size = static_cast<mp_size_t>(sign) * sz;
        mp_limb_t* data = (sz > mpz->_mp_alloc) ? reinterpret_cast<mp_limb_t*>(_mpz_realloc(mpz, sz)) : mpz->_mp_d;
        ar >> boost::serialization::make_array(data, sz);
        // Strangely, mp_size_t is a long int, but the struct member that corresponds is an int.
        // This makes no sense to me.
        mpz->_mp_size = static_cast<int>(signed_size);
#endif
    }
}

/*
 * The GMP float type is a little different. It has four fields: precision (which gives num_limbs = precision + 1);
 * size, which is the number of limbs that are populated; exponent, which is the base of the limbs; and a pointer to
 * limbs. This means our helper functions for mpz_t won't work.
 */

template<typename Archive>
inline void save_to_archive(Archive& ar, mpf_srcptr mpf)
{
    // For some reason, mpf_t doesn't have a lot of helper functions so let's just reach into the struct itself
    int size = mpf->_mp_size;
    int sign = 0;
    if (size > 0) {
        sign = 1;
    }
    else if (size < 0) {
        size = -size;
        sign = -1;
    }

    ar << boost::serialization::make_nvp("precision", mpf->_mp_prec);
    ar << boost::serialization::make_nvp("size", size);
    ar << boost::serialization::make_nvp("sign", sign);
    ar << boost::serialization::make_nvp("exponent", mpf->_mp_exp);
    ar << boost::serialization::make_nvp("mantissa", boost::serialization::make_array(mpf->_mp_d, size));
}

template<typename Archive>
inline void load_from_archive(Archive& ar, mpf_ptr mpf)
{
    // Assume the mpf is initialised
    int size = 0, precision = 0, sign = 0;

    ar >> boost::serialization::make_nvp("precision", precision);
    ar >> boost::serialization::make_nvp("size", size);
    ar >> boost::serialization::make_nvp("sign", sign);

    mpf_set_prec(mpf, precision);// allocates storage for the limbs;

    mpf->_mp_size = sign * size;
    ar >> boost::serialization::make_nvp("exponent", mpf->_mp_exp);
    ar >> boost::serialization::make_nvp("mantissa", boost::serialization::make_array(mpf->_mp_d, size));
}

}



namespace boost {
namespace serialization {

/*
 * OK, so mpz_t is fundamentally an C array of __mpz_structs. So in serialization, the logic is something like as
 * follows:
 *  - Wrap mpz_t into a boost::serialization::array_wrapper<__mpz_struct>
 *  - call serialize method on this wrapper
 *  - dispatch to serialize for __mpz_struct (which doesn't exist)
 *
 *  Unfortunately, this means that implementing serialize for mpz_t does nothing, since the dispatch
 *  is called for boost::serialization::array_wrapper<__mpz_struct>. The solution is to specialize this
 *  wrapper class and provide our own serialization here. This is horrible, I really would prefer to not
 *  do this.
 */

template<>
struct array_wrapper<const __mpz_struct> : public wrapper_traits<const array_wrapper<__mpz_struct>> {
public:
    array_wrapper(const __mpz_struct* t, std::size_t sz) : m_t(t), m_element_count(sz)
    {}

    template<typename Archive>
    void serialize(Archive& ar, const unsigned long int version)
    {
        boost::serialization::split_member(ar, *this, version);
    }

    template<typename Archive>
    void load(Archive& ar, const unsigned long int)
    {
        // Should this really be implemented?
    }

    template<typename Archive>
    void save(Archive& ar, const unsigned long int) const
    {
        gmp_ser_helpers::save_to_archive(ar, m_t);
    }

private:
    const __mpz_struct* const m_t;
    const std::size_t m_element_count;
};

template<>
struct array_wrapper<__mpz_struct> : public wrapper_traits<const array_wrapper<__mpz_struct>> {
public:
    array_wrapper(__mpz_struct* t, std::size_t sz) : m_t(t), m_element_count(sz)
    {}

    template<typename Archive>
    void serialize(Archive& ar, const unsigned long int version)
    {
        boost::serialization::split_member(ar, *this, version);
    }

    template<typename Archive>
    void load(Archive& ar, const unsigned long int)
    {
        mpz_init(m_t);
        gmp_ser_helpers::load_from_archive(ar, m_t);
    }

    template<typename Archive>
    void save(Archive& ar, const unsigned long int) const
    {
        gmp_ser_helpers::save_to_archive(ar, m_t);
    }

private:
    __mpz_struct* const m_t;
    const std::size_t m_element_count;
};

/*
 * We have to use the same trick as above to serialize mpq_t, since it is also an array of __mpq_struct.
 */
template<>
struct array_wrapper<const __mpq_struct> : public wrapper_traits<const array_wrapper<__mpz_struct>> {
public:
    array_wrapper(const __mpq_struct* t, std::size_t sz)
        : m_t(t), m_element_count(sz)
    {}

    template<typename Archive>
    void serialize(Archive& ar, const unsigned long int version)
    {
        boost::serialization::split_member(ar, *this, version);
    }

    template<typename Archive>
    void load(Archive& ar, const unsigned long int)
    {
        // Should this really be implemented?
    }

    template<typename Archive>
    void save(Archive& ar, const unsigned long int) const
    {
        // Save numerator first
        gmp_ser_helpers::save_to_archive(ar, &m_t->_mp_num);
        gmp_ser_helpers::save_to_archive(ar, &m_t->_mp_den);
    }

private:
    const __mpq_struct* const m_t;
    const std::size_t m_element_count;
};

template<>
struct array_wrapper<__mpq_struct> : public wrapper_traits<const array_wrapper<__mpz_struct>> {
public:
    array_wrapper(__mpq_struct* t, std::size_t sz)
        : m_t(t), m_element_count(sz)
    {}

    template<typename Archive>
    void serialize(Archive& ar, const unsigned long int version)
    {
        boost::serialization::split_member(ar, *this, version);
    }

    template<typename Archive>
    void load(Archive& ar, const unsigned long int)
    {
        mpq_init(m_t);
        // Load numerator first
        gmp_ser_helpers::load_from_archive(ar, &m_t->_mp_num);
        gmp_ser_helpers::load_from_archive(ar, &m_t->_mp_den);
    }

    template<typename Archive>
    void save(Archive& ar, const unsigned long int) const
    {
        // Save numerator first
        gmp_ser_helpers::save_to_archive(ar, &m_t->_mp_num);
        gmp_ser_helpers::save_to_archive(ar, &m_t->_mp_den);
    }

private:
    __mpq_struct* const m_t;
    const std::size_t m_element_count;
};

/*
 * Same trick again for mpf_t
 */

template<>
struct array_wrapper<const __mpf_struct> : public wrapper_traits<const array_wrapper<__mpz_struct>> {
public:
    array_wrapper(const __mpf_struct* t, std::size_t sz)
        : m_t(t), m_element_count(sz)
    {}

    template<typename Archive>
    void serialize(Archive& ar, const unsigned long int version)
    {
        boost::serialization::split_member(ar, *this, version);
    }

    template<typename Archive>
    void load(Archive& ar, const unsigned long int)
    {
        // Should this really be implemented?
    }

    template<typename Archive>
    void save(Archive& ar, const unsigned long int) const
    {
        gmp_ser_helpers::save_to_archive(ar, m_t);
    }

private:
    const __mpf_struct* const m_t;
    const std::size_t m_element_count;
};

template<>
struct array_wrapper<__mpf_struct> : public wrapper_traits<const array_wrapper<__mpz_struct>> {
public:
    array_wrapper(__mpf_struct* t, std::size_t sz)
        : m_t(t), m_element_count(sz)
    {}

    template<typename Archive>
    void serialize(Archive& ar, const unsigned long int version)
    {
        boost::serialization::split_member(ar, *this, version);
    }

    template<typename Archive>
    void load(Archive& ar, const unsigned long int)
    {
        mpf_init(m_t);
        gmp_ser_helpers::load_from_archive(ar, m_t);
    }

    template<typename Archive>
    void save(Archive& ar, const unsigned long int) const
    {
        gmp_ser_helpers::save_to_archive(ar, m_t);
    }

private:
    __mpf_struct* const m_t;
    const std::size_t m_element_count;
};

/// Boost::multiprecision GMP backend serialization methods

/*
 * gmp_int is a very simple wrapper around mpz_t.
 */

template<typename Archive>
void save(Archive& ar, const multiprecision::backends::gmp_int& integer_obj, const unsigned int /*version*/)
{
    const mpz_t& mpz = integer_obj.data();
    gmp_ser_helpers::save_to_archive(ar, mpz);
}

template<typename Archive>
void load(Archive& ar, multiprecision::backends::gmp_int& integer_obj, const unsigned int /*version*/)
{
    gmp_ser_helpers::load_from_archive(ar, integer_obj.data());
}

/*
 * gmp_rational is a wrapper around mpq_t, which is an array of 1 element of the struct __mpq_struct containing
 * two __mpz_structs so we can use the same helpers to simply write out the numerator and denominator sequentially
 */

template<typename Archive>
void save(Archive& ar, const multiprecision::backends::gmp_rational& rational_obj, const unsigned int /*version*/)
{
    const mpq_t& mpq = rational_obj.data();
    // Save numerator first
    gmp_ser_helpers::save_to_archive(ar, &mpq->_mp_num);
    gmp_ser_helpers::save_to_archive(ar, &mpq->_mp_den);
}

template<typename Archive>
void load(Archive& ar, multiprecision::backends::gmp_rational& rational_obj, const unsigned int /*version*/)
{
    mpq_t& mpq = rational_obj.data();
    // Load numerator first
    gmp_ser_helpers::load_from_archive(ar, &mpq->_mp_num);
    gmp_ser_helpers::load_from_archive(ar, &mpq->_mp_den);
}

/*
 * gmp_float is a wrapper around an mpf_t. Let's provide an implementation for this too
 */
template<typename Archive, unsigned Digits10>
void serialize(Archive& ar, multiprecision::backends::gmp_float<Digits10>& float_obj, const unsigned int version)
{
    split_free(ar, float_obj, version);
}

template<typename Archive, unsigned Digits10>
void save(Archive& ar, const multiprecision::backends::gmp_float<Digits10>& float_obj, const unsigned int /*version*/)
{
    gmp_ser_helpers::save_to_archive(ar, float_obj.data());
}

template<typename Archive, unsigned Digits10>
void load(Archive& ar, multiprecision::backends::gmp_float<Digits10>& float_obj, const unsigned int /*version*/)
{
    gmp_ser_helpers::load_from_archive(ar, float_obj.data());
}

}} // namespaces


#endif//GMP_SER_GMP_SER_H
