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
 * Sample from a Chi-squared distribution
 */

#ifndef _statslib_rchisq_HPP
#define _statslib_rchisq_HPP

//
// scalar output

template<typename T>
statslib_inline
return_t<T>
rchisq(const T dof_par, rand_engine_t& engine);

template<typename T>
statslib_inline
return_t<T>
rchisq(const T dof_par, const ullint_t seed_val = std::random_device{}());

//
// vector/matrix output

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename mT, typename T1>
statslib_inline
mT
rchisq(const ullint_t n, const ullint_t k, const T1 dof_par);
#endif

//
// include implementation files

#include "rchisq.ipp"

#endif
