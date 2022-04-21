/* *************************************************************

Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurk� and Arend Janssen.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */

//  lie.h

// Include once wrapper
#ifndef DJC_COROPA_LIBALGEBRA_LIEH_SEEN
#define DJC_COROPA_LIBALGEBRA_LIEH_SEEN

namespace alg {

template<typename Coeff>
class lie_multiplication
{

    typedef typename Coeff::SCA scalar_t;

    /// The recursive key product.
    // Caution Reutenauer only states the induction where the state variable
    // is the whole expression (the tree or the pre lie sum of products...)
    // we do not update that state and start over each time.
    template<typename Lie>
    static Lie _prod(typename Lie::KEY const& k1, typename Lie::KEY const& k2)
    {
        typedef Lie lie_t;
        typedef typename Lie::BASIS basis_t;
        typedef typename basis_t::PARENT parent_t;
        typedef typename basis_t::KEY lie_key_t;

        typename Lie::BASIS& basis = Lie::basis;

#ifdef DEBUG// NO LOOPS IN RECURSION FOR _PROD
        static std::map<PARENT, unsigned> counter;
        PARENT tmp(k1, k2);
        if (++counter[tmp] > 1) {
            assert(false);
        };// add debug code here};
        assert(k1 < k2);
#endif// DEBUG

        // We look up for the desired product in our basis.
        parent_t parents(k1, k2);
        typename std::map<parent_t, lie_key_t>::const_iterator it;
        it = lie_t::basis.find(parents);
        if (it != basis.reverse_map_end()) {
            // [k1,k2] exists in the basis.
            return lie_t(it->second);
        }
        else
        // [k1,k2] does not exists in the basis.
        {
            // Since k1 <= k2, k2 is not a letter because if it was a letter,
            // then also k1, which is impossible since [k1,k2] is not in the basis.
            // similarly setting k2 = [k3,k4], it cannot be the case that k3 < k1
            // or this would be a hall tree. So k1 < k3 < k4
            // We use Jacobi: [k1,k2] = [k1,[k3,k4]]] = [[k1,k3],k4]-[[k1,k4],k3]
            lie_key_t k3(basis.lparent(k2));
            lie_key_t k4(basis.rparent(k2));
            lie_t result(prod<lie_t>(k1, k3) * (lie_t)k4);
            result.sub_mul(prod<lie_t>(k1, k4), (lie_t)k3);
            return result;
        }
    }

    template<typename Lie>
    static std::map<typename Lie::BASIS::PARENT, Lie> prime_prod_cache_table()
    {
        typedef typename Lie::BASIS::PARENT parent_t;
        std::map<parent_t, Lie> rv;
        rv[parent_t(0, 0)] = Lie();
        return rv;
    }

    template<typename Lie>
    static const algebra<typename Lie::BASIS, Coeff, lie_multiplication, vectors::sparse_vector>&
    prod(typename Lie::KEY const& k1, typename Lie::KEY const& k2)
    {
        typedef algebra<typename Lie::BASIS, Coeff, lie_multiplication, vectors::sparse_vector> lie_t;

        typedef typename Lie::BASIS basis_t;
        typedef typename basis_t::PARENT parent_t;
        //typedef typename basis_t::KEY lie_key_t;

        static const lie_t zero;
        DEG target_degree = lie_t::basis.degree(k1) + lie_t::basis.degree(k2);// degrees[k1] + degrees[k2];
        if (target_degree > lie_t::BASIS::degree_tag::max_degree) {
            return zero;
        }// degree truncation

        static boost::recursive_mutex table_access;
        static std::map<parent_t, lie_t> table(prime_prod_cache_table<lie_t>());
        // get exclusive recursive access for the thread
        boost::lock_guard<boost::recursive_mutex> lock(table_access);
        // [A,A] = 0.
        if (k1 == k2) {
            return table[parent_t(0, 0)];
        }

        typename std::map<parent_t, lie_t>::iterator it;
        parent_t p(k1, k2);
        it = table.find(p);
        if (it == table.end()) {
            lie_t* ptr = &(table[p] = ((p.first < p.second) ? _prod<lie_t>(k1, k2) : -_prod<lie_t>(k2, k1)));
            return *ptr;
        }
        else {
            return it->second;
        }
    }

    template<typename Lie, typename Transform>
    class index_operator
    {
        Transform m_transform;

        typedef Lie lie_t;
        typedef typename lie_t::const_iterator lie_iterator;

    public:
        void operator()(scalar_t* result_ptr, scalar_t const* lhs_ptr, scalar_t const* rhs_ptr, DIMN const lhs_target,
                        DIMN const rhs_target, bool assign = false)
        {
            scalar_t lhs;
            if (assign) {
                for (IDIMN i = 0; i < static_cast<IDIMN>(lhs_target); ++i) {
                    lhs = lhs_ptr[i];
                    for (IDIMN j = 0; j < static_cast<IDIMN>(rhs_target); ++j) {
                        *(result_ptr++) = m_transform(Coeff::mul(lhs, rhs_ptr[j]));
                    }
                }
            }
            else {
                for (IDIMN i = 0; i < static_cast<IDIMN>(lhs_target); ++i) {
                    lhs = lhs_ptr[i];
                    for (IDIMN j = 0; j < static_cast<IDIMN>(rhs_target); ++j) {
                        *(result_ptr++) += m_transform(Coeff::mul(lhs, rhs_ptr[j]));
                    }
                }
            }
        }
    };

    template<typename Lie, typename Transform>
    class key_operator
    {
        Transform m_transform;

    public:
#if __cplusplus >= 201103UL

        template<typename... Args>
        explicit key_operator(Args&&... args)
            : m_transform(std::forward<Args>(args)...)
        {}

#else
        template<typename Arg>
        explicit key_operator(Arg arg) : m_transform(arg)
        {}
#endif

        template<typename Vector>
        void
        operator()(Vector& result, typename Vector::KEY const& lhs_key, scalar_t const& lhs_val,
                   typename Vector::KEY const& rhs_key, scalar_t const& rhs_val)
        {
            result.add_scal_prod(prod<Lie>(lhs_key, rhs_key), m_transform(Coeff::mul(lhs_val, rhs_val)));
        }
    };

public:
    template<typename Algebra, typename Operator>
    Algebra& multiply_and_add(Algebra& result, Algebra const& lhs, Algebra const& rhs, Operator op) const
    {
        key_operator<Algebra, Operator> kt(op);
        lhs.buffered_apply_binary_transform(result, rhs, kt);
        return result;
    }

    template<typename Algebra, typename Operator>
    Algebra&
    multiply_and_add(Algebra& result, Algebra const& lhs, Algebra const& rhs, Operator op, DEG const max_depth) const
    {
        key_operator<Algebra, Operator> kt(op);
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
        key_operator<Algebra, Operator> kt(op);
        lhs.unbuffered_apply_binary_transform(rhs, kt);
        return lhs;
    }

    template<typename Algebra, typename Operator>
    Algebra& multiply_inplace(Algebra& lhs, Algebra const& rhs, Operator op, DEG const max_depth) const
    {
        key_operator<Algebra, Operator> kt(op);
        lhs.unbuffered_apply_binary_transform(rhs, kt, max_depth);
        return lhs;
    }
};

/**
 * @brief  A specialisation of the algebra class with a Lie basis.
 *
 * Mathematically, the algebra of Lie instances is a free Lie associative
 * algebra. With respect to the inherited algebra class, the essential
 * distinguishing feature of this class is the basis class used, and in
 * particular the basis::prod() member function. Thus, the most important
 * information is in the definition of lie_basis. Notice that this associative
 * algebra of lie elements does not includes as a sub-algebra the associative
 * algebra corresponding to the SCALAR type. In other words, only the scalar
 * zero corresponds to a Lie element (the zero one) which is the neutral
 * element of the addition operation. There is no neutral element for the
 * product (free Lie product).
 */
template<typename Coeff, DEG n_letters, DEG max_degree,
         template<typename, typename, typename...> class VectorType,
         typename... Args>
class lie : public algebra<
                    lie_basis<n_letters, max_degree>, Coeff, lie_multiplication<Coeff>, VectorType, Args...>
{
    typedef lie_multiplication<Coeff> multiplication_t;

public:
    /// The basis type.
    typedef lie_basis<n_letters, max_degree> BASIS;
    /// Import of the KEY type.
    typedef typename BASIS::KEY KEY;
    /// The algebra type.
    typedef algebra<BASIS, Coeff, multiplication_t, VectorType, Args...> ALG;
    /// The sparse_vector type.
    typedef typename ALG::VECT VECT;

    typedef typename Coeff::SCA SCA;
    typedef typename Coeff::RAT RAT;

    /// Import of the iterator type.
    typedef typename ALG::iterator iterator;
    /// Import of the constant iterator type.
    typedef typename ALG::const_iterator const_iterator;

public:
    /// Default constructor. Zero lie element.
    lie()
    {}

    /// Copy constructor.
    lie(const lie& l)
        : ALG(l)
    {}

    /// Constructs an instance from an algebra instance.
    lie(const ALG& a)
        : ALG(a)
    {}

    /// Constructs an instance from a sparse_vector instance.
    lie(const VECT& v)
        : ALG(v)
    {}

    /// Constructs a unidimensional instance from a given key (with scalar one).
    explicit lie(const KEY& k)
        : ALG(k)
    {}

    // explicit lie(LET letter, const SCA& s)
    //	// flawed as basis is possibly not yet constructed
    //	: ALG(VECT::basis.keyofletter(letter), s) {}
    /// Constructs a unidimensional instance from a key and a scalar.
    explicit lie(const KEY& k, const SCA& s)
        : ALG(k, s)
    {}

    /// Construct an instance from pointers to data range
    lie(SCA const* begin,
        SCA const* end) : ALG(begin, end)
    {
    }

    lie& operator=(const lie&) = default;
    //lie& operator=(lie&&) noexcept = default;

public:
    /// Replaces the occurrences of letters in s by Lie elements in v.
    inline friend lie replace(const lie& src, const std::vector<LET>& s, const std::vector<const lie*>& v)
    {
        lie result;
        std::map<KEY, lie> table;
        const_iterator i;
        for (i = src.begin(); i != src.end(); ++i) {
            result.add_scal_prod(VECT::basis.replace(i->key(), s, v, table), i->value());
        }
        return result;
    }


#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
private:

    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, unsigned int const /* version */) {
        ar & boost::serialization::base_object<ALG>(*this);
    }
#endif
};

}// namespace alg

// Include once wrapper
#endif// DJC_COROPA_LIBALGEBRA_LIEH_SEEN

// EOF.
