# third_party/random

## What are these files?

They are copies of `random` and `algorithm` headers from the [llvm](https://github.com/llvm-mirror/libcxx/tree/master/include) standard library.


## Motivation

Users complained about stability of random numbers across machines when setting seed across different platforms. See issue [#379](https://github.com/grf-labs/grf/issues/379).

As [pointed out](https://github.com/grf-labs/grf/issues/379#issuecomment-480641123) by @jtibshirani:

> the mersenne twister has the same implementation across platforms, the other random methods may differ from compiler to compiler


In PR [#469](https://github.com/grf-labs/grf/pull/469), @halflearned included this reduced copy of the relevant headers, ensuring random number generation is done in a consistent way across compilers.


## How to reproduce

#### File `random.hpp`

Extract only the relevant classes and functions from `random.hpp`:

+ `__independent_bits_engine`
+ `__clz`
+ `__log2`
+ `__log2_imp`
+ `__generate_canonical`
+ `uniform_int_distribution`
+ `poisson_distribution`
+ `uniform_real_distribution`
+ `normal_distribution`
+ `exponential_distribution`
+ `discrete_distribution`

From each class, remove all methods associated with `operator<<`.

Find and remove the following `_LIBCPP` macros:
+ `_LIBCPP_BEGIN_NAMESPACE_STD`
+ `_LIBCPP_CONSTEXPR`
+ `_LIBCPP_DISABLE_UBSAN_UNSIGNED_INTEGER_CHECK`
+ `_LIBCPP_END_NAMESPACE_STD`
+ `_LIBCPP_HAS_NO_PRAGMA_SYSTEM_HEADER`
+ `_LIBCPP_INLINE_VISIBILITY`
+ `_LIBCPP_MSVCRT`
+ `_LIBCPP_POP_MACROS`
+ `_LIBCPP_PUSH_MACROS`
+ `_LIBCPP_RANDOM`
+ `_LIBCPP_TEMPLATE_VIS`
+ `_LIBCPP_TYPE_VIS`
+ `_LIBCPP_USING_DEV_RANDOM`

Find and replace prefix:
+ `_VSTD::` -> `std::`

Add `namespace nonstd`.


#### File `algorithm.hpp`

Include modified `random.hpp`

Extract relevant class:

+ `shuffle`

Replace prefix:

+ `std::uniform_int_distribution` -> `nonstd::uniform_int_distribution`

Add `namespace nonstd`.
