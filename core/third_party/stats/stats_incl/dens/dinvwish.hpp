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
 * pdf of the inverse-Wishart distribution
 */

#ifndef _statslib_dinvwish_HPP
#define _statslib_dinvwish_HPP

#ifdef STATS_ENABLE_MATRIX_FEATURES

template<typename mT, typename pT, typename not_arma_mat<mT>::type* = nullptr>
statslib_inline
return_t<pT> dinvwish(const mT& X, const mT& Psi_par, const pT nu_par, const bool log_form = false);

// specializations
#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename pT>
statslib_inline
eT dinvwish(const ArmaMat<eT>& X, const ArmaMat<eT>& Psi_par, const pT nu_par, const bool log_form = false);
#endif

#include "dinvwish.ipp"

#endif

#endif
