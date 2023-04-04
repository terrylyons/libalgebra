/* *************************************************************

Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurko and Arend Janssen, Sam Morley

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */

//  area_tensor_multiplication.h

// Include once wrapper
#ifndef DJC_COROPA_LIBALGEBRA_area_TENSOR_MULTIPLICATIONH_SEEN
#define DJC_COROPA_LIBALGEBRA_area_TENSOR_MULTIPLICATIONH_SEEN


#include "tensor.h"

#if 0
#include "half_shuffle_tensor_basis.h"
#include "area_tensor_basis.h"

namespace alg {

template<typename Coeff>
class area_tensor_multiplication
{

    typedef typename Coeff::SCA scalar_t;

    /// Computes recursively the half shuffle product of two keys
    template<typename Tensor>
    static Tensor _prod(typename Tensor::KEY const& k1, typename Tensor::KEY const& k2)
    {
        typedef typename Tensor::BASIS basis_t;
        typedef half_shuffle_tensor<Coeff, basis_t::s_no_letters, basis_t::s_max_degree> half_shuffle_tensor_t;

        Tensor result;

        if ((basis_t::s_max_degree == 0) || (k1.size() + k2.size() <= basis_t::s_max_degree)) {
            // is this undefined behaviour - probably it is defined but risky behaviour
            // the conversion is from an r value so there can be no question of coexistence of pointers and aliasing
            // The standard says the compiler will not check that the conversion makes sense and that one should determine that it is of this type by some other means.

            // Commented this out, it is equivalent to below.
            //const half_shuffle_tensor_t& temp = (half_shuffle_tensor_t(k2) * half_shuffle_tensor_t(k1)
            //	- half_shuffle_tensor_t(k1) * half_shuffle_tensor_t(k2));
            //typename half_shuffle_tensor_t::const_iterator i;
            //for (i = temp.begin(); i != temp.end(); ++i) {
            //	result[i->key()] += i->value();
            //}

            static_assert(sizeof(Tensor) == sizeof(half_shuffle_tensor_t), "Tensor and half_shuffle_tensor_t must be binary compatible");
            auto tmp = half_shuffle_tensor_t(k2) * half_shuffle_tensor_t(k1)
                    - half_shuffle_tensor_t(k1) * half_shuffle_tensor_t(k2);

            for (auto& t : tmp) {
                result.add_scal_prod(t.key(), t.value());
            }
        }
        return result;
    }

    /// The half shuffle product of two basis elements.
    /**
    Returns the area_tensor obtained by 
                  k1 \prec k2 = k1.lparent()*(k1.rparent() \prec k2 + k2 \prec k1.rparent()) 
    where * is the normal concatenation or tensor product. 
    If k1 is zero it is zero, and if k2 is zero and k1 is not zero it is k1.
    The already computed products are stored in a static multiplication 
    table to speed up further calculations.

    If it were not for the caching of multiplication one might remove the recursion
    */
    template<typename Tensor>
    static const Tensor& prod(typename Tensor::KEY const& k1, typename Tensor::KEY const& k2)
    {
        typedef typename Tensor::KEY key_t;
        static boost::recursive_mutex table_access;
        // get exclusive recursive access for the thread
        boost::lock_guard<boost::recursive_mutex> lock(table_access);

        typedef std::map<std::pair<key_t, key_t>, Tensor> TABLE_T;
        static TABLE_T table;
        typename TABLE_T::iterator it;
        std::pair<key_t, key_t> p(k1, k2);
        it = table.find(p);
        if (it == table.end()) {
            return table[p] = _prod<Tensor>(k1, k2);
        }
        else {
            return it->second;
        }
    }

    template<typename Transform>
    class index_operator
    {
        Transform m_transform;

    public:
        index_operator(Transform t) : m_transform(t) {}

        void operator()(scalar_t* result_ptr, scalar_t const* lhs_ptr, scalar_t const* rhs_ptr, DIMN const lhs_target,
                        DIMN const rhs_target, bool assign = false) {}
    };

    template<typename Transform>
    class key_operator
    {
        Transform m_transform;

    public:
        key_operator(Transform t) : m_transform(t) {}

        template<typename Vector>
        void
        operator()(Vector& result, typename Vector::KEY const& lhs_key, scalar_t const& lhs_val,
                   typename Vector::KEY const& rhs_key, scalar_t const& rhs_val)
        {
            result.add_scal_prod(prod<Vector>(lhs_key, rhs_key), m_transform(Coeff::mul(lhs_val, rhs_val)));
        }
    };

public:
    template<typename Algebra, typename Operator>
    Algebra& multiply_and_add(Algebra& result, Algebra const& lhs, Algebra const& rhs, Operator op) const
    {
        key_operator<Operator> kt(op);
        lhs.buffered_apply_binary_transform(result, rhs, kt);
        return result;
    }

    template<typename Algebra, typename Operator>
    Algebra&
    multiply_and_add(Algebra& result, Algebra const& lhs, Algebra const& rhs, Operator op, DEG const max_depth) const
    {
        key_operator<Operator> kt(op);
        lhs.buffered_apply_binary_transform(result, rhs, kt, max_depth);
        return result;
    }

    template<typename Algebra, typename Operator>
    Algebra multiply(Algebra const& lhs, Algebra const& rhs, Operator op) const
    {
        Algebra result;
        multiply_and_add(result, lhs, rhs, op);
        return result;
    }

    template<typename Algebra, typename Operator>
    Algebra multiply(Algebra const& lhs, Algebra const& rhs, Operator op, DEG const max_depth) const
    {
        Algebra result;
        multiply_and_add(result, lhs, rhs, op, max_depth);
        return result;
    }

    template<typename Algebra, typename Operator>
    Algebra& multiply_inplace(Algebra& lhs, Algebra const& rhs, Operator op) const
    {
        key_operator<Operator> kt(op);
        lhs.unbuffered_apply_binary_transform(rhs, kt);
        return lhs;
    }

    template<typename Algebra, typename Operator>
    Algebra& multiply_inplace(Algebra& lhs, Algebra const& rhs, Operator op, DEG const max_depth) const
    {
        key_operator<Operator> kt(op);
        lhs.unbuffered_apply_binary_transform(rhs, kt, max_depth);
        return lhs;
    }
};

/// A specialisation of the algebra class with a half_shuffletensor basis.
/**
   Mathematically, the algebra of area_tensor instances is a shuffle
   associative algebra associated to a free associative algebra. With respect
   to the inherited algebra class, the essential distinguishing feature of
   this class is the basis class used, and in particular the basis::prod()
   member function. Thus, the most important information is in the definition
   of area_tensor_basis. Notice that this associative algebra of free
   tensors includes as a sub-algebra the associative algebra corresponding to
   the SCALAR type. This is permitted by the existence of empty keys in
   area_tensor_basis.
 */
template<typename Coeff, DEG n_letters, DEG max_degree>
class area_tensor : public alg::algebra<
                            area_tensor_basis<n_letters, max_degree>, Coeff, area_tensor_multiplication<Coeff>

                            >
{
    typedef area_tensor_multiplication<Coeff> multiplication_t;

public:
    /// The basis type.
    typedef area_tensor_basis<n_letters, max_degree> BASIS;
    /// Import of the KEY type.
    typedef typename BASIS::KEY KEY;
    /// The algebra type.
    typedef algebra<BASIS, Coeff, multiplication_t> ALG;

    typedef typename Coeff::SCA SCA;
    typedef typename Coeff::RAT RAT;

    /// The sparse_vector type.
    typedef typename ALG::VECT VECT;

    /// Import of the iterator type.
    typedef typename ALG::iterator iterator;
    /// Import of the constant iterator type.
    typedef typename ALG::const_iterator const_iterator;

public:
    /// Default constructor.
    area_tensor(void) {}

    /// Copy constructor.
    area_tensor(const area_tensor& t) : ALG(t) {}

    /// Constructs an instance from a free_tensor instance.
    area_tensor(const free_tensor<Coeff, n_letters, max_degree>& t)
    {
        typename free_tensor<Coeff, n_letters, max_degree>::const_iterator i;
        for (i = t.begin(); i != t.end(); ++i) {
            (*this)[i->key()] += i->value();
        }
    }

    /// Constructs an instance from an algebra instance.
    area_tensor(const ALG& a) : ALG(a) {}

    /// Constructs an instance from a sparse_vector instance.
    area_tensor(const VECT& v) : ALG(v) {}

    /// Constructs a unidimensional instance from a letter and a scalar.
    area_tensor(LET
                        letter,
                const SCA& s)
        : ALG(VECT::basis
                      .keyofletter(letter),
              s)
    {
    }

    /// Constructs a unidimensional instance from a key (basis element).
    explicit area_tensor(const KEY& k) : ALG(k) {}

    /// Constructs a unidimensional instance from a scalar.
    explicit area_tensor(const SCA& s) : ALG(VECT::basis.empty_key, s) {}

public:
    //
    ///// Ensures that the return type is a area_tensor.
    //inline __DECLARE_BINARY_OPERATOR (area_tensor, *, *=, SCA)
    //
    ///// Ensures that the return type is a area_tensor.
    //inline __DECLARE_BINARY_OPERATOR (area_tensor, /, /=, RAT)
    //
    ///// Ensures that the return type is a area_tensor.
    //inline __DECLARE_BINARY_OPERATOR (area_tensor, *, *=, area_tensor)
    //
    ///// Ensures that the return type is a area_tensor.
    //inline __DECLARE_BINARY_OPERATOR (area_tensor, +, +=, area_tensor)
    //
    ///// Ensures that the return type is a area_tensor.
    //inline __DECLARE_BINARY_OPERATOR (area_tensor, -, -=, area_tensor)
    //
    ///// Ensures that the return type is a area_tensor.
    //inline __DECLARE_UNARY_OPERATOR (area_tensor, -, -, ALG)
};

} // namespace alg
#endif
// Include once wrapper
#endif // DJC_COROPA_LIBALGEBRA_area_TENSOR_MULTIPLICATIONH_SEEN

// EOF.
