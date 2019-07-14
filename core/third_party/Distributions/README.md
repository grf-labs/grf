# Distributions

These functions were adapted from `Boost`, version 1.70.0 (released April 12th, 2019 06:04 GMT).

+ Function `uniform_draw` was adapted from `boost/random/uniform_real_distribution.hpp:L33-47`. 
+ Function `poisson_draw` was adapted from `boost/random/poisson_distribution.hpp:L252-332`.

Most of the code is simple copy-and-paste and the logic was kept the same as in the original `boost` code. However, 
some variables were renamed for readability, and some free variables (e.g. `[poisson_table](https://github.com/grf-labs/grf/blob/f6179f72d560404d8dc95fc304eda1274f4eaff5/core/third_party/Distributions/Distributions.h#L64)` 
are now defined within the function.) 