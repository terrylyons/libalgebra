#ifndef half_shuffle_tensor_basis_h__
#define half_shuffle_tensor_basis_h__

#include "implementation_types.h"

namespace alg {

/// The monoid of words of a finite number of letters with the half-shuffle product.
/**
This is the basis used to implement the shuffle_tensor class as a
specialization of the algebra class. This basis is the Free Associative
Algebra basis with a finite number of letters, with the shuffle product.
The shuffle_tensor_basis is a container of keys. A key is the
implementation of a word of letters. The prod() member function
corresponds to the half-shuffle product of the two keys given as arguments.
This product is neither associative or commutative. Letters can be seen as
particular basis keys, i.e. words of length one. The empty word is a
special key used for the embedding of letters (words of length one).
*/
template<DEG n_letters, DEG max_degree>
class half_shuffle_tensor_basis : public alg::tensor_basis<n_letters, max_degree>,
                                  public alg::base_basis<With_Degree, n_letters,
                                                           max_degree>
{
public:
    /// The tensor_basis type.
    typedef tensor_basis<n_letters, max_degree> TBASIS;
    /// Import of the KEY type.
    typedef typename TBASIS::KEY KEY;
    /// Import of the MAP type.
    // typedef typename TBASIS::MAP MAP;

    typedef alg::basis::with_degree<max_degree> degree_tag;
    typedef alg::basis::ordered<std::less<KEY>> ordering_tag;

public:
    /// Default constructor.
    half_shuffle_tensor_basis(void) {}
};
//#endif

namespace vectors {

template<DEG n_letters, DEG max_depth, typename Field>
struct vector_type_selector<half_shuffle_tensor_basis<n_letters, max_depth>, Field> {
    typedef half_shuffle_tensor_basis<n_letters, max_depth> BASIS;
    typedef typename BASIS::KEY KEY;
    typedef sparse_vector<BASIS, Field,
#ifndef ORDEREDMAP
                          MY_UNORDERED_MAP<KEY, typename Field::S, typename KEY::hash>
#else
                           std::map<KEY, typename Field::S>
#endif
                          >
            type;
};

}// namespace vectors

} // namespace alg

#endif // half_shuffle_tensor_basis_h__
