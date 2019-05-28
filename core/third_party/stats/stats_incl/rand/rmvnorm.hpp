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
 * Sample from a multivariate normal distribution
 */

#ifndef _statslib_rmvnorm_HPP
#define _statslib_rmvnorm_HPP

#ifdef STATS_ENABLE_MATRIX_FEATURES

// vector draw
template<typename vT, typename mT, typename not_arma_mat<mT>::type* = nullptr>
statslib_inline
vT rmvnorm(const vT& mu_par, const mT& Sigma_par, const bool pre_chol = false);

// n samples
template<typename vT, typename mT, typename not_arma_mat<mT>::type* = nullptr>
statslib_inline
mT rmvnorm(const ullint_t n, const vT& mu_par, const mT& Sigma_par, const bool pre_chol = false);

// specializations
#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename mT, typename eT>
statslib_inline
mT rmvnorm(const mT& mu_par, const ArmaMat<eT>& Sigma_par, const bool pre_chol = false);

template<typename mT, typename eT>
statslib_inline
ArmaMat<eT> rmvnorm(const ullint_t n, const mT& mu_par, const ArmaMat<eT>& Sigma_par, const bool pre_chol = false);
#endif

#include "rmvnorm.ipp"

#endif

#endif
