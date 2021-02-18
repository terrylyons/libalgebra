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

/// A specialisation of the algebra class with a free tensor basis.
/**
   Mathematically, the algebra of free_tensor instances is a free associative
   algebra. With respect to the inherited algebra class, the essential
   distinguishing feature of this class is the basis class used, and in
   particular the basis::prod() member function. Thus, the most important
   information is in the definition of free_tensor_basis. Notice that this
   associative algebra of free tensors includes as a sub-algebra the
   associative algebra corresponding to the SCALAR type. This is permitted by
   the existence of empty keys in free_tensor_basis.
 */
template<typename SCA, typename RAT, DEG n_letters, DEG max_degree>
class free_tensor : public algebra<free_tensor_basis<SCA, RAT, n_letters, max_degree> >
{
public:
	/// The basis type.
	typedef free_tensor_basis<SCA, RAT, n_letters, max_degree> BASIS;
	/// Import of the KEY type.
	typedef typename BASIS::KEY KEY;
    /// The algebra type.
    typedef algebra<BASIS> ALG;
	/// The sparse_vector type.
	typedef typename ALG::VECT VECT;

	/// Import of the iterator type.
	typedef typename ALG::iterator iterator;
	/// Import of the constant iterator type.
	typedef typename ALG::const_iterator const_iterator;
public:
	/// Default constructor.
	free_tensor(void) {}
	/// Copy constructor.
	free_tensor(const free_tensor& t)
		: ALG(t) {}
	/// Constructs an instance from a shuffle_tensor instance.
	free_tensor(const shuffle_tensor<SCA, RAT, n_letters, max_degree>& t)
	{
		typename shuffle_tensor<SCA, RAT, n_letters, max_degree>::const_iterator i;
		for (i = t.begin(); i != t.end(); ++i)
			(*this)[i->first] += i->second;
	}
	/// Constructs an instance from an algebra instance.
	free_tensor(const ALG& a)
		: ALG(a) {}
	/// Constructs an instance from a sparse_vector instance.
	free_tensor(const VECT& v)
		: ALG(v) {}	
	/// Constructs a unidimensional instance from a letter and a scalar.
	free_tensor(LET letter, const SCA& s)
		: ALG(VECT::basis.keyofletter(letter), s) {}
	/// Explicit unidimensional constructor from a given key (basis element).
	explicit free_tensor(const KEY& k)
		: ALG(k) {}
	/// Explicit unidimensional constructor from a given scalar.
	explicit free_tensor(const SCA& s)
		: ALG(VECT::basis.empty_key, s) {}
public:
	/// Ensures that the return type is a free_tensor.
  inline __DECLARE_BINARY_OPERATOR(free_tensor,*,*=,SCA)
	/// Ensures that the return type is a free_tensor.
  inline __DECLARE_BINARY_OPERATOR(free_tensor,/,/=,RAT)
	/// Ensures that the return type is a free_tensor.
  inline __DECLARE_BINARY_OPERATOR(free_tensor,*,*=,free_tensor)
	/// Ensures that the return type is a free_tensor.
  inline __DECLARE_BINARY_OPERATOR(free_tensor,+,+=,free_tensor)
	/// Ensures that the return type is a free_tensor.
  inline __DECLARE_BINARY_OPERATOR(free_tensor,-,-=,free_tensor)
	/// Ensures that the return type is a free_tensor.
  inline __DECLARE_UNARY_OPERATOR(free_tensor,-,-,ALG)
	/// Computes the truncated exponential of a free_tensor instance.
	inline friend free_tensor exp(const free_tensor& arg)
	{
		// Computes the truncated exponential of arg
		// 1 + arg + arg^2/2! + ... + arg^n/n! where n = max_degree
		KEY kunit;
		free_tensor result(kunit);
		for (DEG i = max_degree; i >= 1; --i)
		{
			result.mul_scal_div(arg, (RAT)i);
			result += (free_tensor)kunit;
		}
		return result;
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
		if (it != x.end())
			x.erase(it);
		free_tensor result;
		for (DEG i = max_degree; i >= 1; --i)
		{
			if (i % 2 == 0)
				result.sub_scal_div(tunit, (RAT)i);
			else
				result.add_scal_div(tunit, (RAT)i);
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
		// = a^(-1) - x/a*[a^(-1)(1 - x/a + x^2/a^2 + ... + (-1)^(n) x^(n-1)/a^(n-1)))].
		// S_n = a^(-1) + z S_{n-1}; z = - x/a ; S_0 = a^(-1)
		// max_degree must be > 0

		static KEY kunit;
		SCA a(0);
		free_tensor x, z(a);

		const_iterator it(arg.find(kunit));
		if (it == arg.end())
			// const term a is 0;
			throw "divide-by-zero";
		else
		{
			a = (*it).value();
			x = arg;
			x.erase(kunit);
		}

		//S_n = a + z S_{ n - 1 }; z = -x / a; S_0 = a
		//
		// the nonzero scalar component a of the tensor arg restored to a tensor
		free_tensor free_tensor_a_inverse(SCA(1)/a), result(free_tensor_a_inverse);
		// z := - x/a
		z.sub_scal_div(x, a);
		// the iteration
		for (DEG i = 0; i != max_degree; ++i)
			result = free_tensor_a_inverse + z * result;
		return result;
	}
	/// Computes the truncated inverse of a free_tensor instance.
	inline friend free_tensor reflect(const free_tensor& arg)
	{
		// Computes the alternating reflection of arg up to degree max_degree
		// For group-like elements this is the same as the inverse 
		free_tensor ans(SCA(0));
		for (const_iterator it = arg.begin(); it != arg.end(); ++it)
		{
			KEY old_key = it->key();
			SCA old_value = it->value();
			ans[old_key.reverse()] = (old_key.size()%2) ? SCA(0) - old_value : old_value;
		}
		return ans;
	}
};

/// A specialisation of the algebra class with a shuffle tensor basis.
/**
   Mathematically, the algebra of shuffle_tensor instances is a shuffle
   associative algebra associated to a free associative algebra. With respect
   to the inherited algebra class, the essential distinguishing feature of
   this class is the basis class used, and in particular the basis::prod()
   member function. Thus, the most important information is in the definition
   of shuffle_tensor_basis. Notice that this associative algebra of free
   tensors includes as a sub-algebra the associative algebra corresponding to
   the SCALAR type. This is permitted by the existence of empty keys in
   shuffle_tensor_basis.
 */
template<typename SCA, typename RAT, DEG n_letters, DEG max_degree>
class shuffle_tensor : public algebra<shuffle_tensor_basis<SCA, RAT, n_letters, max_degree> >
{
public:
	/// The basis type.
	typedef shuffle_tensor_basis<SCA, RAT, n_letters, max_degree> BASIS;
	/// Import of the KEY type.
	typedef typename BASIS::KEY KEY;
    /// The algebra type.
    typedef algebra<BASIS> ALG;

	/// The sparse_vector type.
	typedef typename ALG::VECT VECT;

	/// Import of the iterator type.
	typedef typename ALG::iterator iterator;
	/// Import of the constant iterator type.
	typedef typename ALG::const_iterator const_iterator;
public:
	/// Default constructor.
	shuffle_tensor(void) {}
	/// Copy constructor.
	shuffle_tensor(const shuffle_tensor& t)
		: ALG(t) {}
	/// Constructs an instance from a free_tensor instance.
	shuffle_tensor(const free_tensor<SCA, RAT, n_letters, max_degree>& t)
	{
		typename free_tensor<SCA, RAT, n_letters, max_degree>::const_iterator i;
		for (i = t.begin(); i != t.end(); ++i)
			(*this)[i->first] += i->second;
	}
	/// Constructs an instance from an algebra instance.
	shuffle_tensor(const ALG& a)
		: ALG(a) {}
	/// Constructs an instance from a sparse_vector instance.
	shuffle_tensor(const VECT& v)
		: ALG(v) {}	
	/// Constructs a unidimensional instance from a letter and a scalar.
	shuffle_tensor(LET letter, const SCA& s)	
		: ALG(VECT::basis.keyofletter(letter), s) {}
	/// Constructs a unidimensional instance from a key (basis element).
	explicit shuffle_tensor(const KEY& k)
		: ALG(k) {}
	/// Constructs a unidimensional instance from a scalar.
	explicit shuffle_tensor(const SCA& s)	
		: ALG(VECT::basis.empty_key, s) {}
public:
	/// Ensures that the return type is a shuffle_tensor.
  inline __DECLARE_BINARY_OPERATOR(shuffle_tensor,*,*=,SCA)
	/// Ensures that the return type is a shuffle_tensor.
  inline __DECLARE_BINARY_OPERATOR(shuffle_tensor,/,/=,RAT)
	/// Ensures that the return type is a shuffle_tensor.
  inline __DECLARE_BINARY_OPERATOR(shuffle_tensor,*,*=,shuffle_tensor)
	/// Ensures that the return type is a shuffle_tensor.
  inline __DECLARE_BINARY_OPERATOR(shuffle_tensor,+,+=,shuffle_tensor)
	/// Ensures that the return type is a shuffle_tensor.
  inline __DECLARE_BINARY_OPERATOR(shuffle_tensor,-,-=,shuffle_tensor)
	/// Ensures that the return type is a shuffle_tensor.
  inline __DECLARE_UNARY_OPERATOR(shuffle_tensor,-,-,ALG)
};

// Include once wrapper
#endif // DJC_COROPA_LIBALGEBRA_TENSORH_SEEN

//EOF.
