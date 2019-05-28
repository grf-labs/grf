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
 * Sample from a logistic distribution
 */

#ifndef _statslib_rlogis_HPP
#define _statslib_rlogis_HPP

//
// scalar output

template<typename T1, typename T2>
statslib_inline
common_return_t<T1,T2>
rlogis(const T1 mu_par, const T2 sigma_par, rand_engine_t& engine);

template<typename T1, typename T2>
statslib_inline
common_return_t<T1,T2>
rlogis(const T1 mu_par, const T2 sigma_par, const ullint_t seed_val = std::random_device{}());

//
// vector/matrix output

#ifdef STATS_ENABLE_MATRIX_FEATURES
template<typename mT, typename T1, typename T2>
statslib_inline
mT
rlogis(const ullint_t n, const ullint_t k, const T1 mu_par, const T2 sigma_par);
#endif

//
// include implementation files

#include "rlogis.ipp"

#endif
