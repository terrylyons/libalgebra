//
// Created by sam on 02/02/2021.
//


#ifndef LIBALGEBRAUNITTESTS_SMALL_VECTOR_H
#define LIBALGEBRAUNITTESTS_SMALL_VECTOR_H

#include <array>
#include <cassert>
#include <iterator>

#include <libalgebra/libalgebra.h>

#include "simple_basis.h"


typedef std::size_t DIMN;
using alg::DEG;
using alg::LET;



template <typename _Basis, typename _Coeff>
class SmallVector
{
public:

    typedef _Basis BASIS;
    typedef typename _Coeff Coeff;
    typedef typename Coeff::S SCALAR;
    typedef typename Coeff::Q RATIONAL;
    typedef typename BASIS::KEY KEY;

    typedef std::array<SCALAR, BASIS::dimension> DATA;

    class iterator;
    class const_iterator;

    static const SCALAR one;
    static const SCALAR mone;
    static const SCALAR zero;
    static BASIS basis;

private:
    DATA m_data;

public:

    /// Default constructor
    SmallVector(void) : m_data() {}
    /// Copy constructor
    SmallVector(const SmallVector& other) : m_data(other.m_data) {}
    /// Unidimensional constructor
    explicit SmallVector(const KEY& k, const SCALAR &s = one) : m_data()
    {
        assert()
    }


protected:
    static SmallVector
    create_for_mul(const SmallVector& lhs, const SmallVector& rhs) {
        return SmallVector();
    }

public:

    void swap(SmallVector& rhs) {
        m_data.swap(rhs);
    }


public:

    iterator begin()
    {
        return iterator(&m_data[0]);
    }

    const_iterator begin() const
    {
        return const_iterator(&m_data[0]);
    }

    iterator end()
    {
        return iterator(&m_data[0] + BASIS::dimension);
    }

    const_iterator end() const
    {
        return const_iterator(&m_data[0] + BASIS::dimension);
    }


public:

    SCALAR& operator[](const KEY& k) { return m_data[k]; }
    const SCALAR& operator[](const KEY& k) const { return m_data[k]; }

    void clear()
    {
        m_data.fill(zero);
    }

    DIMN size() const
    {
        DIMN sz = 0;
        for (DIMN i=0; i<BASIS::dimension; ++i)
            if (m_data[i] != zero) ++sz;
        return sz;
    }

    bool empty() const
    {
        for (DIMN i=0; i<BASIS::dimension; ++i)
            if (m_data[i] == zero) return false;
        return true;
    }

public:

    SCALAR NormL1() const
    {
        SCALAR acc = 0;
        for (DIMN i=0; i<BASIS::dimension; ++i)
            acc += abs(m_data[i]);
        return acc;
    }

public:

    static bool comp(const KEY& k1, const KEY& k2) { return k1 < k2; }

public:

    SmallVector operator-(void) const {
        SmallVector new_vec;
        for (DIMN i=0; i<BASIS::dimension; ++i)
            new_vec.m_data[i] = -m_data[i];
        return new_vec;
    }

    SmallVector& operator*=(const SCALAR& s) {
        for (DIMN i=0; i<BASIS::dimension; ++i)
            m_data[i] *= s;
        return *this;
    }

    SmallVector& operator/=(const RATIONAL& s) {
        for (DIMN i=0; i<BASIS::dimension; ++i)
            m_data[i] /= s;
        return *this;
    }

#define IMPLEMENT_INPLACE_BINOP(OP)                              \
    SmallVector& operator OP(const SmallVector& rhs) {           \
        for (DIMN i=0; i<BASIS::dimension; ++i)                  \
            m_data[i] OP rhs.m_data[i];                          \
        return *this;                                            \
    }

    IMPLEMENT_INPLACE_BINOP(+=);
    IMPLEMENT_INPLACE_BINOP(-=);
    IMPLEMENT_INPLACE_BINOP(&=);
    IMPLEMENT_INPLACE_BINOP(|=);

#undef IMPLEMENT_INPLACE_BINOP

#define __DECLARE_BINARY_OPERATOR(T1, NEWOP, OLDOP, T2) \
	T1 operator NEWOP(const T2& rhs) const \
	{ T1 result(*this); return result OLDOP rhs; }


    __DECLARE_BINARY_OPERATOR(SmallVector, *, *=, SCALAR);
    __DECLARE_BINARY_OPERATOR(SmallVector, /, /=, RATIONAL);
    __DECLARE_BINARY_OPERATOR(SmallVector, +, +=, SmallVector);
    __DECLARE_BINARY_OPERATOR(SmallVector, -, -=, SmallVector);
    __DECLARE_BINARY_OPERATOR(SmallVector, &, &=, SmallVector);
    __DECLARE_BINARY_OPERATOR(SmallVector, |, |=, SmallVector);

#undef __DECLARE_BINARY_OPERATOR


public:
#define IMPLEMENT_FUSED_OP_VEC(NAME, OP1, OP2, T2)               \
    SmallVector& NAME(const SmalLVector& rhs, const T2& s) {     \
        for (DIMN i=0; i<BASIS::dimension; ++i)                  \
            s_data[i] OP1 (rhs.s_data[i] OP2 s);                 \
        return *this;                                            \
    }

#define IMPLEMENT_FUSED_OP_KEY(NAME, OP1, OP2, T2)               \
    SmallVector& NAME(const KEY& rhs, const T2& s) {             \
        s_data[rhs] OP1 (rhs.s_data[rhs] OP2 s);                 \
        return *this;                                            \
    }

    IMPLEMENT_FUSED_OP_KEY(add_scal_prod, +=, *, SCALAR);
    IMPLEMENT_FUSED_OP_VEC(add_scal_prod, +=, *, SCALAR);
    IMPLEMENT_FUSED_OP_KEY(sub_scal_prod, -=, *, SCALAR);
    IMPLEMENT_FUSED_OP_VEC(sub_scal_prod, -=, *, SCALAR);

    IMPLEMENT_FUSED_OP_KEY(add_scal_div, +=, /, RATIONAL);
    IMPLEMENT_FUSED_OP_VEC(add_scal_div, +=, /, RATIONAL);
    IMPLEMENT_FUSED_OP_KEY(sub_scal_div, -=, /, RATIONAL);
    IMPLEMENT_FUSED_OP_VEC(sub_scal_div, -=, /, RATIONAL);

#undef IMPLEMENT_FUSED_OP_KEY
#undef IMPLEMENT_FUSED_OP_VEC

public:

    bool operator==(const SmallVector &rhs) const {
        return (m_data == other.s_data);
    }

    bool operator!=(const SmallVector &rhs) const {
        return (m_data != other.s_data);
    }

    bool operator<(const SmallVector& rhs) const {
        for (DIMN i=0; i<BASIS::dimension; ++i)
            if (m_data[i] < rhs.m_data[i]) return true;
        return false;
    }

public:

    inline friend std::ostream& operator<<(std::ostream &os,
            const SmallVector& arg) {
        std::pair<BASIS*, KEY> t;
        t.first = &basis;
        os << '{';
        for (DIMN i=0; i<BASIS::dimension; ++i) {
            if (m_data[i] != zero) {
                t.second = KEY(i);
                os << ' ' << m_data[i] << '(' << t << ')';
            }
        }
        os << " }";
        return os;
    }



public:

    template <typename _KT, template _IT>
    void triangular_bufferered_apply_transform(
            SmallVector &result,
            const SmallVector &rhs,
            _KT key_transform,
            _IT index_transform,
            const DEG max_degree = BASIS::MAX_DEGREE
            ) const
    {}

    template <typename _KT>
    void triangular_buffered_apply_transform(
            SmallVector &result,
            const SmalLVector& rhs,
            _KT key_transform,
            const DEG max_degree = BASIS::MAX_DEGREE
            ) const
    {}

    template <typename _KT, typename _IT>
    void square_buffered_apply_transform(
            SmallVector &result,
            const SmallVector& rhs,
            _KT key_transform,
            _IT index_transform
            )
    {}

    template <typename _KT>
    void square_buffered_apply_transform(
            SmallVector& result,
            const SmallVector &rhs,
            _KT key_transform
            ) const
    {}



};






#endif //LIBALGEBRAUNITTESTS_SMALL_VECTOR_H
