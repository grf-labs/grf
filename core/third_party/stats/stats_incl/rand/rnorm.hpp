/*################################################################################
  ##
  ##   Copyright (C) 2011-2019 Keith O'Hara
  ##
  ##   This file is part of the StatsLib C++ library.
  ##
  ##   Licensed under the Apache License, Version 2.0 (the "License");
  ##   you may not use this file except in compliance with the License.
  ##   You may obtain a copy of the License at
  ##
  ##       http://www.apache.org/licenses/LICENSE-2.0
  ##
  ##   Unless required by applicable law or agreed to in writing, software
  ##   distributed under the License is distributed on an "AS IS" BASIS,
  ##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  ##   See the License for the specific language governing permissions and
  ##   limitations under the License.
  ##
  ################################################################################*/

/*
 * Sample from a normal distribution
 */

#ifndef _statslib_rnorm_HPP
#define _statslib_rnorm_HPP

//
// scalar output

template<typename T1, typename T2>
statslib_inline
common_return_t<T1,T2>
rnorm(const T1 mu_par, const T2 sigma_par, rand_engine_t& engine);

template<typename T1, typename T2>
statslib_inline
common_return_t<T1,T2>
rnorm(const T1 mu_par, const T2 sigma_par, const ullint_t seed_val = std::random_device{}());

template<typename T = double>
statslib_inline
T rnorm();

//
// vector/matrix output

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename mT, typename T1 = double, typename T2 = double>
statslib_inline
mT
rnorm(const ullint_t n, const ullint_t k, const T1 mu_par = T1(0), const T2 sigma_par = T2(1));
#endif

//
// include implementation files

#include "rnorm.ipp"

#endif
