# Libalgebra
Libalgebra is a header only template library for performing high level computations in 
abstract mathematical structures, such as the tensor algebra and free Lie algebras. The core
components in the library are the vector classes. In the current version of libalgebra,
there are three vector templates: `dense_vector`, `sparse_vector`, and `hybrid_vector`.
These underlying vector types are unified in a single vector interface class `vector`, 
upon which the rest of the class hierarchy is based.

A object of the vector template class describes a generic element in the linear span
of a basis, provided by a template argument, with coefficients from a field, also
provided by a template argument. For example, an element of the tensor algebra with 4
letters with real coefficients (described by double precision floating point numbers),
would be represented by objects of the class
```c++
alg::vector<alg::tensor_basis<4, 2>, alg::double_field>
```
The second template argument to `tensor_basis` describes the maximum degree of tensor words
to consider in the basis.

Built on top of the vector interface is the general algebra template class. The 
algebra class takes an additional template argument in addition to those required
by `vector` that describes the multiplication operation. This is further specialised
to give top level template classes for the free tensor algebras, shuffle tensor
algebras, free Lie algebras, polynomial algebras, and poly-Lie algebras. Also included
is a class that provides an interface for using the Campbell-Baker-Hausdorff formula
for computing the solution to an equation of the form exp(x) = exp(a)exp(b) in the
free Lie algebras, and a class for converting between various algebra types.

## Installation
This library uses CMake to manage builds. There are several ways to include Libalgebra
into your projects. The simplest method is to include libalgebra as a git submodule, or
even just clone into a subdirectory of your main project folder. A better way is to
use the CMake `FetchContent` module to obtain and make available the libalgebra target
to your own CMake project file. To do this, you can use the commands
```cmake
include(FetchContent)

FetchContent_Declare(libalgebra
        GIT_REPOSITORY https://github.com/terrylyons/libalgebra.git
        GIT_TAG master
        )
FetchContent_MakeAvailable(libalgebra)
```
Once this is done, you can link the `Libalgebra::Libalgebra` target in your own targets.

Libalgebra has dependencies on Boost and GMP (or MPIR on Windows), so these will need to
be installed and available for CMake to find; we use the `FindBoost` and `FindBigNum` CMake
modules to locate and link the required components.
