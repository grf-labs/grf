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
 * cdf of the Poisson distribution
 */

#ifndef _statslib_ppois_HPP
#define _statslib_ppois_HPP

//
// single input

template<typename T>
statslib_constexpr
return_t<T> ppois(const llint_t x, const T rate_par, const bool log_form = false) noexcept;

//
// vector/matrix input

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename rT = common_return_t<eT,T1>>
statslib_inline
std::vector<rT>
ppois(const std::vector<eT>& x, const T1 rate_par, const bool log_form = false);
#endif

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename rT = common_return_t<eT,T1>>
statslib_inline
ArmaMat<rT>
ppois(const ArmaMat<eT>& X, const T1 rate_par, const bool log_form = false);

template<typename mT, typename tT, typename T1>
statslib_inline
mT 
ppois(const ArmaGen<mT,tT>& X, const T1 rate_par, const bool log_form = false);
#endif

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename rT = common_return_t<eT,T1>, bool To = blaze::columnMajor>
statslib_inline
BlazeMat<rT,To>
ppois(const BlazeMat<eT,To>& X, const T1 rate_par, const bool log_form = false);
#endif

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename rT = common_return_t<eT,T1>, int iTr = Eigen::Dynamic, int iTc = Eigen::Dynamic>
statslib_inline
EigenMat<rT,iTr,iTc>
ppois(const EigenMat<eT,iTr,iTc>& X, const T1 rate_par, const bool log_form = false);
#endif

//
// include implementation files

#include "ppois.ipp"

#endif
