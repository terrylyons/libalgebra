//
// Created by sam on 12/02/2021.
//


/*
 * Multiplication operator for simple_integer_basis
 */
struct pointwise_multiplication
{
    typedef unsigned key_type;

    template <typename Transform>
    struct key_operator
    {
        key_operator(Transform fn) : m_transform(fn)
        {}

        template <typename Vector, typename Scalar>
        void operator()(
                Vector& result, const key_type& k1, const Scalar& s1, const key_type& k2, const Scalar& s2)
        {
            if (k1 == k2) {
                result.add_scal_prod(k1, m_transform(s1 * s2));
            }
        }

    private:
        Transform m_transform;
    };

    template <typename Algebra, typename Operator>
    Algebra &multiply_and_add(Algebra &result, Algebra const &lhs, Algebra const &rhs, Operator op) const
    {
        key_operator <Operator> kt(op);
        lhs.buffered_apply_binary_transform(result, rhs, kt);
        return result;
    }

    template <typename Algebra, typename Operator> Algebra &
    multiply_and_add(Algebra &result, Algebra const &lhs, Algebra const &rhs, Operator op, DEG const max_depth) const
    {
        key_operator <Operator> kt(op);
        lhs.buffered_apply_binary_transform(result, rhs, kt, max_depth);
        return result;
    }

    template <typename Algebra, typename Operator>
    Algebra multiply(Algebra const &lhs, Algebra const &rhs, Operator op) const
    {
        Algebra result;
        multiply_and_add(result, lhs, rhs, op);
        return result;
    }

    template <typename Algebra, typename Operator>
    Algebra multiply(Algebra const &lhs, Algebra const &rhs, Operator op, DEG const max_depth) const
    {
        Algebra result;
        multiply_and_add(result, lhs, rhs, op, max_depth);
        return result;
    }

    template <typename Algebra, typename Operator>
    Algebra &multiply_inplace(Algebra &lhs, Algebra const &rhs, Operator op) const
    {
        key_operator <Operator> kt(op);
        lhs.unbuffered_apply_binary_transform(rhs, kt);
        return lhs;
    }

    template <typename Algebra, typename Operator>
    Algebra &multiply_inplace(Algebra &lhs, Algebra const &rhs, Operator op, DEG const max_depth) const
    {
        key_operator <Operator> kt(op);
        lhs.unbuffered_apply_binary_transform(rhs, kt, max_depth);
        return lhs;
    }
};


