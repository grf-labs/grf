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
 * quantile function of the Weibull distribution
 */

#ifndef _statslib_qweibull_HPP
#define _statslib_qweibull_HPP

//
// single input

template<typename T1, typename T2, typename T3>
statslib_constexpr
common_return_t<T1,T2,T3>
qweibull(const T1 p, const T2 shape_par, const T3 scale_par) noexcept;

//
// vector/matrix input

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT = common_return_t<eT,T1,T2>>
statslib_inline
std::vector<rT>
qweibull(const std::vector<eT>& x, const T1 shape_par, const T2 scale_par);
#endif

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT = common_return_t<eT,T1,T2>>
statslib_inline
ArmaMat<rT>
qweibull(const ArmaMat<eT>& X, const T1 shape_par, const T2 scale_par);

template<typename mT, typename tT, typename T1, typename T2>
statslib_inline
mT
qweibull(const ArmaGen<mT,tT>& X, const T1 shape_par, const T2 scale_par);
#endif

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT = common_return_t<eT,T1,T2>, bool To = blaze::columnMajor>
statslib_inline
BlazeMat<rT,To>
qweibull(const BlazeMat<eT,To>& X, const T1 shape_par, const T2 scale_par);
#endif

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT = common_return_t<eT,T1,T2>, int iTr = Eigen::Dynamic, int iTc = Eigen::Dynamic>
statslib_inline
EigenMat<rT,iTr,iTc>
qweibull(const EigenMat<eT,iTr,iTc>& X, const T1 shape_par, const T2 scale_par);
#endif

//
// include implementation files

#include "qweibull.ipp"

#endif
