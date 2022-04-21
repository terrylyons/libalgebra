/* *************************************************************

Copyright 2010 Terry Lyons, Stephen Buckley, Djalil Chafai,
Greg Gyurk� and Arend Janssen.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */

//  tensor.h

// Include once wrapper
#ifndef DJC_COROPA_LIBALGEBRA_TENSORH_SEEN
#define DJC_COROPA_LIBALGEBRA_TENSORH_SEEN

namespace alg {

template<typename Coeff>
class free_tensor_multiplication;

template<typename Coeff>
class shuffle_tensor_multiplication;

template<typename Coeff, DEG n_letters, DEG max_degree, typename...>
class shuffle_tensor;

template<typename Coeff>
class free_tensor_multiplication
{

    typedef typename Coeff::SCA scalar_t;

    template<typename Transform>
    class index_operator
    {
        Transform m_transform;

    public:
        index_operator(Transform t)
            : m_transform(t)
        {}

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

    template<typename Transform>
    class key_operator
    {
        Transform m_transform;

    public:
        key_operator(Transform t)
            : m_transform(t)
        {}

        template<typename Vector>
        void
        operator()(Vector& result, typename Vector::KEY const& lhs_key, scalar_t const& lhs_val,
                   typename Vector::KEY const& rhs_key, scalar_t const& rhs_val)
        {
            result.add_scal_prod(lhs_key * rhs_key, m_transform(Coeff::mul(lhs_val, rhs_val)));
        }
    };

public:
    template<typename Algebra, typename Operator>
    Algebra& multiply_and_add(Algebra& result, Algebra const& lhs, Algebra const& rhs, Operator op) const
    {
        key_operator<Operator> kt(op);
        index_operator<Operator> it(op);
        lhs.buffered_apply_binary_transform(result, rhs, kt, it);
        return result;
    }

    template<typename Algebra, typename Operator>
    Algebra&
    multiply_and_add(Algebra& result, Algebra const& lhs, Algebra const& rhs, Operator op, DEG const max_depth) const
    {
        key_operator<Operator> kt(op);
        index_operator<Operator> it(op);
        lhs.buffered_apply_binary_transform(result, rhs, kt, it, max_depth);
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
        index_operator<Operator> it(op);
        lhs.unbuffered_apply_binary_transform(rhs, kt, it);
        return lhs;
    }

    template<typename Algebra, typename Operator>
    Algebra& multiply_inplace(Algebra& lhs, Algebra const& rhs, Operator op, DEG const max_depth) const
    {
        key_operator<Operator> kt(op);
        index_operator<Operator> it(op);
        lhs.unbuffered_apply_binary_transform(rhs, kt, it, max_depth);
        return lhs;
    }
};

/**
 * @brief A specialisation of the algebra class with a free tensor basis.
 *
 * Mathematically, the algebra of free_tensor instances is a free associative
 * algebra. With respect to the inherited algebra class, the essential
 * distinguishing feature of this class is the basis class used, and in
 * particular the basis::prod() member function. Thus, the most important
 * information is in the definition of free_tensor_basis. Notice that this
 * associative algebra of free tensors includes as a sub-algebra the
 * associative algebra corresponding to the SCALAR type. This is permitted by
 * the existence of empty keys in free_tensor_basis.
 */
template<typename Coeff, DEG n_letters, DEG max_degree,
         template<typename, typename, typename...> class VectorType,
         typename... Args>
class free_tensor : public algebra<
                            free_tensor_basis<n_letters, max_degree>, Coeff, free_tensor_multiplication<Coeff>, VectorType, Args...>
{
    typedef free_tensor_multiplication<Coeff> multiplication_t;

public:
    /// The basis type.
    typedef free_tensor_basis<n_letters, max_degree> BASIS;
    /// Import of the KEY type.
    typedef typename BASIS::KEY KEY;
    /// The algebra type.
    typedef algebra<BASIS, Coeff, multiplication_t, VectorType, Args...> ALG;

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
    free_tensor()
    {}

    /// Copy constructor.
    free_tensor(const free_tensor& t)
        : ALG(t)
    {}

    /// Constructs an instance from a shuffle_tensor instance.
    free_tensor(const shuffle_tensor<Coeff, n_letters, max_degree>& t)
    {
        typename shuffle_tensor<Coeff, n_letters, max_degree>::const_iterator i;
        for (i = t.begin(); i != t.end(); ++i) {
            (*this)[i->key()] += i->value();
        }
    }

    /// Constructs an instance from an algebra instance.
    free_tensor(const ALG& a)
        : ALG(a)
    {}

    /// Constructs an instance from a sparse_vector instance.
    free_tensor(const VECT& v)
        : ALG(v)
    {}

    /// Constructs a unidimensional instance from a letter and a scalar.
    free_tensor(LET letter, const SCA& s) : ALG(VECT::basis.keyofletter(letter), s)
    {
    }

    /// Constructor from pointers to data range
    free_tensor(SCA const* begin,
                SCA const* end) : ALG(begin, end)
    {
    }

    /// Explicit unidimensional constructor from a given key (basis element).
    explicit free_tensor(const KEY& k)
        : ALG(k)
    {}

    /// Explicit unidimensional constructor from a given scalar.
    explicit free_tensor(const SCA& s)
        : ALG(VECT::basis.empty_key, s)
    {}

    /// Constructor from pointers to data range with offset
    free_tensor(DIMN
                        offset,
                SCA const* begin, SCA const* end) : ALG(offset, begin, end)
    {
    }

    free_tensor& operator=(const free_tensor&) = default;
    //free_tensor& operator=(free_tensor&&) noexcept = default;

public:
    /// Ensures that the return type is a free_tensor.
    inline free_tensor operator*(const SCA& rhs) const
    {
        free_tensor result(*this);
        result *= rhs;
        return result;
    }

    /// Ensures that the return type is a free_tensor.
    inline free_tensor operator/(const RAT& rhs) const
    {
        free_tensor result(*this);
        result /= rhs;
        return result;
    }

    /// Ensures that the return type is a free_tensor.
    inline free_tensor operator*(const free_tensor& rhs) const
    {
        free_tensor result(*this);
        result *= rhs;
        return result;
    }

    /// Ensures that the return type is a free_tensor.
    inline free_tensor operator+(const free_tensor& rhs) const
    {
        free_tensor result(*this);
        result += rhs;
        return result;
    }

    /// Ensures that the return type is a free_tensor.
    inline free_tensor operator-(const free_tensor& rhs) const
    {
        free_tensor result(*this);
        result -= rhs;
        return result;
    }

    /// Ensures that the return type is a free_tensor.
    inline free_tensor operator-() const
    {
        return free_tensor(ALG::operator-());
    }

    /// Computes the truncated exponential of a free_tensor instance.
    inline friend free_tensor exp(const free_tensor& arg)
    {
        // Computes the truncated exponential of arg
        // 1 + arg + arg^2/2! + ... + arg^n/n! where n = max_degree
        KEY kunit;
        free_tensor result(kunit);
        for (DEG i = max_degree; i >= 1; --i) {
            result.mul_scal_div(arg, (RAT)i);
            result += (free_tensor)
                    kunit;
        }
        return result;
    }

    /**
     * Fused multiply exponential operation for free tensors.
     *
     * Computes a*exp(x) using a modified Horner's method for the case when x does
     * not have a constant term. If the argument exp_arg has a constant term, it
     * is ignored.
     *
     * For a real number x, one can expand exp(x) up to degree n as
     *
     *     1 + b_1 x(1 + b_2 x(1 + ... b_n x(1)) ...)
     *
     * where each b_i has the value 1/i. This formulae works when x is a free
     * tensor, or indeed any element in an unital (associative) algebra. Working
     * through the result of multiplying on the left by another element a in the
     * above gives the expansion
     *
     *     a + b1 (a + b_2 (a + ... b_n (a)x) ... x)x.
     *
     * This is the result of a*exp(x). In a non-commutative algebra this need not
     * be equal to exp(x)*a.
     *
     * @param exp_arg free_tensor (const reference) to expentiate (x).
     * @return free_tensor containing a*exp(x)
     */
    free_tensor fmexp(const free_tensor& exp_arg) const
    {
        free_tensor result(*this), x(exp_arg);
        KEY kunit;
        typename free_tensor::iterator unit_elt;

        if ((unit_elt = x.find(kunit)) != x.end() && unit_elt->value() != VECT::zero) {
            x.erase(unit_elt);
        }

        for (DEG i = max_degree; i >= 1; --i) {
            result.mul_scal_div(x, static_cast<RAT>(i), max_degree - i + 1);
            result += *this;
        }

        return result;
    }

    /// Inplace version of fmexp
    free_tensor& fmexp_inplace(const free_tensor& exp_arg)
    {
        free_tensor self(*this), x(exp_arg);
        KEY kunit;
        typename free_tensor::iterator unit_elt;

        if ((unit_elt = x.find(kunit)) != x.end() && unit_elt->value() != VECT::zero) {
            x.erase(unit_elt);
        }

        for (DEG i = max_degree; i >= 1; --i) {
            this->mul_scal_div(x, static_cast<RAT>(i), max_degree - i + 1);
            *this += self;
        }

        return *this;
    }

    /// Computes the truncated logarithm of a free_tensor instance.
    inline friend free_tensor log(const free_tensor& arg)
    {
        // Computes the truncated log of arg up to degree max_degree
        // The coef. of the constant term (empty word in the monoid) of arg
        // is forced to 1.
        // log(arg) = log(1+x) = x - x^2/2 + ... + (-1)^(n+1) x^n/n.
        // max_degree must be > 0
        KEY kunit;
        free_tensor tunit(kunit);
        free_tensor x(arg);
        iterator it = x.find(kunit);
        if (it != x.end()) {
            x.erase(it);
        }
        free_tensor result;

        for (DEG i = max_degree; i >= 1; --i) {
            if (i % 2 == 0) {
                result.sub_scal_div(tunit, (RAT)i);
            }
            else {
                result.add_scal_div(tunit, (RAT)i);
            }
            result *= x;
        }

        return result;
    }

    /// Computes the truncated inverse of a free_tensor instance.
    inline friend free_tensor inverse(const free_tensor& arg)
    {
        // Computes the truncated inverse of arg up to degree max_degree
        // An exception is thrown if the leading term is zero.
        // the module assumes
        // (a+x)^(-1) = (a(1+x/a))^(-1)
        //  = a^(-1)(1 - x/a + x^2/a^2 + ... + (-1)^(n) x^n/a^n)
        // = a^(-1) - x/a*[a^(-1)(1 - x/a + x^2/a^2 + ... + (-1)^(n)
        // x^(n-1)/a^(n-1)))]. S_n = a^(-1) + z S_{n-1}; z = - x/a ; S_0 = a^(-1)
        // max_degree must be > 0

        static KEY kunit;
        SCA a(0);
        free_tensor x, z(a);

        const_iterator it(arg.find(kunit));
        if (it == arg.end()) {
            // const term a is 0;
            throw std::invalid_argument("divide-by-zero");
        }
        else {
            a = (*it).value();
            x = arg;
            x.erase(kunit);
        }

        // S_n = a + z S_{ n - 1 }; z = -x / a; S_0 = a
        //
        //  the nonzero scalar component a of the tensor arg restored to a tensor
        free_tensor free_tensor_a_inverse(SCA(1) / a), result(free_tensor_a_inverse);
        // z := - x/a
        z.sub_scal_div(x, a);
        // the iteration
        for (DEG i = 0; i != max_degree; ++i) {
            result = free_tensor_a_inverse + z * result;
        }
        return result;
    }

    /// Computes the truncated inverse of a free_tensor instance.
    inline friend free_tensor reflect(const free_tensor& arg)
    {
        // Computes the alternating reflection of arg up to degree max_degree
        // For group-like elements this is the same as the inverse
        free_tensor ans(SCA(0));
        for (const_iterator it = arg.begin(); it != arg.end(); ++it) {
            KEY old_key = it->key();
            SCA old_value = it->value();
            ans[old_key.reverse()] = (old_key.size() % 2) ? SCA(0) - old_value : old_value;
        }
        return ans;
    }

#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
private:
    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive& ar, unsigned int const /* version */)
    {
        ar& boost::serialization::base_object<ALG>(*this);
    }
#endif
};

template<typename Coeff>
class shuffle_tensor_multiplication
{

    typedef typename Coeff::SCA scalar_t;

    /// Computes recursively the shuffle product of two keys
    template<typename Tensor>
    static Tensor _prod(typename Tensor::KEY const& k1, typename Tensor::KEY const& k2)
    {
        typedef typename Tensor::KEY key_type;

        typedef typename Tensor::BASIS basis_t;

        Tensor result;
        // unsigned i, j;
        const scalar_t one(+1);

        if ((basis_t::degree_tag::max_degree == 0) || (k1.size() + k2.size() <= basis_t::degree_tag::max_degree)) {
            if (k1.size() == 0) {
                result[k2] = one;
                return result;
            }
            if (k2.size() == 0) {
                result[k1] = one;
                return result;
            }
            // k1.size() >= 1 and k2.size() >= 1
            // let's just implement the multiplication
            const Tensor& first = prod<Tensor>(k1.rparent(), k2);
            const Tensor& second = prod<Tensor>(k1, k2.rparent());
            const key_type k1l{k1.lparent()}, k2l{k2.lparent()};

            typename Tensor::const_iterator cit;

            for (cit = first.begin(); cit != first.end(); ++cit) {
                result[k1l * cit->key()] += Coeff::one;
            }
            for (cit = second.begin(); cit != second.end(); ++cit) {
                result[k2l * cit->key()] += Coeff::one;
            }
        }

        return result;
    }

    /// The shuffle product of two basis elements.
    /**
    Returns the shuffle_tensor obtained by the concatenation product of two
    keys viewed as words of letters. The result is a unidimensional
    shuffle_tensor with a unique key (the concatenation of k1 and k2)
    associated to the +1 scalar. The already computed products are stored in
    a static mutiplication table to speed up further calculations.
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
        std::pair<key_t, key_t> p(std::min(k1, k2), std::max(k1, k2));
        it = table.find(p);
        if (it == table.end()) {
            return table[p] = _prod<Tensor>(k1, k2);
        }
        else {
            return it->second;
        }
    }

    template<typename Transform>
    class key_operator
    {
        Transform m_transform;

    public:
        explicit key_operator(Transform t)
            : m_transform(t)
        {}

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

/**
 * @brief A specialisation of the algebra class with a shuffle tensor basis.
 *
 * Mathematically, the algebra of shuffle_tensor instances is a shuffle
 * associative algebra associated to a free associative algebra. With respect
 * to the inherited algebra class, the essential distinguishing feature of
 * this class is the basis class used, and in particular the basis::prod()
 * member function. Thus, the most important information is in the definition
 * of shuffle_tensor_basis. Notice that this associative algebra of free
 * tensors includes as a sub-algebra the associative algebra corresponding to
 * the SCALAR type. This is permitted by the existence of empty keys in
 * shuffle_tensor_basis.
 */
template<typename Coeff, DEG n_letters, DEG max_degree, typename...>
class shuffle_tensor : public algebra<
                               shuffle_tensor_basis<n_letters, max_degree>, Coeff, shuffle_tensor_multiplication<Coeff>>
{
    typedef shuffle_tensor_multiplication<Coeff> multiplication_t;

public:
    /// The basis type.
    typedef shuffle_tensor_basis<n_letters, max_degree> BASIS;
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
    using ALG::ALG;

    /// Default constructor.
    shuffle_tensor()
    {}

    /// Copy constructor.
    shuffle_tensor(const shuffle_tensor& t)
        : ALG(t)
    {}

    /// Constructs an instance from a free_tensor instance.
    shuffle_tensor(const free_tensor<Coeff, n_letters, max_degree>& t)
    {
        typename free_tensor<Coeff, n_letters, max_degree>::const_iterator i;
        for (i = t.begin(); i != t.end(); ++i) {
            (*this)[i->key()] += i->value();
        }
    }

    /// Constructs an instance from an algebra instance.
    shuffle_tensor(const ALG& a)
        : ALG(a)
    {}

    /// Constructs an instance from a sparse_vector instance.
    shuffle_tensor(const VECT& v)
        : ALG(v)
    {}

    /// Constructs a unidimensional instance from a letter and a scalar.
    shuffle_tensor(LET
                           letter,
                   const SCA& s)
        : ALG(VECT::basis
                      .keyofletter(letter),
              s)
    {
    }

    /// Constructs a unidimensional instance from a key (basis element).
    explicit shuffle_tensor(const KEY& k)
        : ALG(k)
    {}

    /// Constructs a unidimensional instance from a scalar.
    explicit shuffle_tensor(const SCA& s)
        : ALG(VECT::basis.empty_key, s)
    {}

    shuffle_tensor& operator=(const shuffle_tensor&) = default;
    shuffle_tensor& operator=(shuffle_tensor&&) noexcept = default;

public:
    /// Ensures that the return type is a shuffle_tensor.
    inline shuffle_tensor operator*(const SCA& rhs) const
    {
        shuffle_tensor result(*this);
        result *= rhs;
        return result;
    }

    /// Ensures that the return type is a shuffle_tensor.
    inline shuffle_tensor operator/(const RAT& rhs) const
    {
        shuffle_tensor result(*this);
        result /= rhs;
        return result;
    }

    /// Ensures that the return type is a shuffle_tensor.
    inline shuffle_tensor operator*(const shuffle_tensor& rhs) const
    {
        shuffle_tensor result(*this);
        result *= rhs;
        return result;
    }

    /// Ensures that the return type is a shuffle_tensor.
    inline shuffle_tensor operator+(const shuffle_tensor& rhs) const
    {
        shuffle_tensor result(*this);
        result += rhs;
        return result;
    }

    /// Ensures that the return type is a shuffle_tensor.
    inline shuffle_tensor operator-(const shuffle_tensor& rhs) const
    {
        shuffle_tensor result(*this);
        result -= rhs;
        return result;
    }

    /// Ensures that the return type is a shuffle_tensor.
    inline shuffle_tensor operator-(void) const { return shuffle_tensor(ALG::operator-()); }

#ifdef LIBALGEBRA_ENABLE_SERIALIZATION
private:
    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive& ar, unsigned int const /* version */)
    {
        ar& boost::serialization::base_object<ALG>(*this);
    }
#endif
};

}// namespace alg
// Include once wrapper
#endif// DJC_COROPA_LIBALGEBRA_TENSORH_SEEN

// EOF.
