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
 * pdf of the Multivariate Normal distribution
 */

#ifndef _statslib_dmvnorm_HPP
#define _statslib_dmvnorm_HPP

#ifdef STATS_ENABLE_MATRIX_FEATURES
template<typename vT, typename mT, typename eT = double>
statslib_inline
eT dmvnorm(const vT& X, const vT& mu_par, const mT& Sigma_par, bool log_form = false);

#include "dmvnorm.ipp"
#endif

#endif
