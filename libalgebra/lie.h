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



template <typename Lie>
class lie_multiplication {


public:
    typedef Lie lie_t;
    typedef typename lie_t::KEY lie_key_t;


private:
    typedef typename lie_t::BASIS basis_t;
    typedef typename lie_t::PARENT parent_t;

    /// The recursive key product.
    // Caution Reutenauer only states the induction where the state variable
    // is the whole expression (the tree or the pre lie sum of products...)
    // we do not update that state and start over each time.
    lie_t _prod(const lie_key_t &k1, const lie_key_t &k2) {
#ifdef DEBUG // NO LOOPS IN RECURSION FOR _PROD
        static std::map<PARENT, unsigned> counter;
        PARENT tmp(k1, k2);
        if (++counter[tmp] > 1) { assert(false); }; // add debug code here};
        assert(k1 < k2);
#endif // DEBUG

        // We look up for the desired product in our basis.
        parent_t parents(k1, k2);
        typename std::map<parent_t, lie_key_t>::const_iterator it;
        it = lie_t::basis.reverse_map.find(parents);
        if (it != lie_t::basis.reverse_map.end()) {
            // [k1,k2] exists in the basis.
            return lie_t(it->second);
        } else
            // [k1,k2] does not exists in the basis.
        {
            // Since k1 <= k2, k2 is not a letter because if it was a letter,
            // then also k1, which is impossible since [k1,k2] is not in the basis.
            // similarly setting k2 = [k3,k4], it cannot be the case that k3 < k1
            // or this would be a hall tree. So k1 < k3 < k4
            // We use Jacobi: [k1,k2] = [k1,[k3,k4]]] = [[k1,k3],k4]-[[k1,k4],k3]
            lie_key_t k3(lparent(k2));
            lie_key_t k4(rparent(k2));
            lie_t result(prod(k1, k3) * (lie_t) k4);
            result.sub_mul(prod(k1, k4), (lie_t) k3);
            return result;
        }
    }

    static std::map<parent_t, lie_t> prime_prod_cache_table() {
        std::map<parent_t, lie_t> rv;
        rv[parent_t(0, 0)] = lie_t();
        return rv;
    }

    inline const lie_t &prod(const lie_key_t &k1, const lie_key_t &k2) {
        static const lie_t zero;
        DEG target_degree = lie_t::basis.degree(k1) + lie_t::basis.degree(k2); //degrees[k1] + degrees[k2];
        if (target_degree > lie_t::basis.max_degree) {
            return zero;
        } // degree truncation

        // We grow up the basis up to the desired degree.
        lie_t::basis.growup(target_degree);

        static boost::recursive_mutex table_access;
        static std::map<parent_t, lie_t> table(prime_prod_cache_table());
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
            lie_t *ptr = &(table[p] = ((p.first < p.second) ? _prod(k1, k2) : -_prod(k2, k1)));
            return *ptr;
        } else {
            return it->second;
        }
    }


};




/// A specialisation of the algebra class with a Lie basis.
/**
   Mathematically, the algebra of Lie instances is a free Lie associative
   algebra. With respect to the inherited algebra class, the essential
   distinguishing feature of this class is the basis class used, and in
   particular the basis::prod() member function. Thus, the most important
   information is in the definition of lie_basis. Notice that this associative
   algebra of lie elements does not includes as a sub-algebra the associative
   algebra corresponding to the SCALAR type. In other words, only the scalar
   zero corresponds to a Lie element (the zero one) which is the neutral
   element of the addition operation. There is no neutral element for the
   product (free Lie product).
 */
template <typename Coeff, DEG n_letters, DEG max_degree, typename VectorType>
class lie : public algebra<lie_basis < n_letters, max_degree>,
        Coeff,
        lie_multiplication < lie<Coeff, n_letters, max_degree, VectorType> >,
        VectorType> {
public:
/// The basis type.
typedef lie_basis <n_letters, max_degree> BASIS;
/// Import of the KEY type.
typedef typename BASIS::KEY KEY;
/// The algebra type.
typedef algebra <BASIS, Coeff, VectorType> ALG;
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
lie(void) {}

/// Copy constructor.
lie(const lie &l) : ALG(l) {}

/// Constructs an instance from an algebra instance.
lie(const ALG &a) : ALG(a) {}

/// Constructs an instance from a sparse_vector instance.
lie(const VECT &v) : ALG(v) {}

/// Constructs a unidimensional instance from a given key (with scalar one).
explicit lie(const KEY &k) : ALG(k) {}

///// Constructs a unidimensional instance from a letter and a scalar.
//explicit lie(LET letter, const SCA& s)
//	// flawed as basis is possibly not yet constructed
//	: ALG(VECT::basis.keyofletter(letter), s) {}
/// Constructs a unidimensional instance from a key and a scalar.
explicit lie(const KEY &k, const SCA &s) : ALG(k, s) {}

public:

/// Replaces the occurrences of letters in s by Lie elements in v.
inline friend lie replace(const lie &src, const std::vector<LET> &s, const std::vector<const lie *> &v) {
    lie result;
    std::map<KEY, lie> table;
    const_iterator i;
    for (i = src.begin(); i != src.end(); ++i) {
        result.add_scal_prod(VECT::basis.replace(i->key(), s, v, table), i->value());
    }
    return result;
}
};

// Include once wrapper
#endif // DJC_COROPA_LIBALGEBRA_LIEH_SEEN

//EOF.
