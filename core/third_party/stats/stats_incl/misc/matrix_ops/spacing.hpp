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
 * for internal use only; used to switch between the different matrix libraries
 */

//
// get the spacing dimensions

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT>
statslib_inline
ullint_t
spacing(const std::vector<eT>& X)
{
    STATS_UNUSED_PAR(X);
    return ullint_t(1);
}
#endif

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT>
statslib_inline
ullint_t
spacing(const ArmaMat<eT>& X)
{
    return n_cols(X);
}
#endif

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, bool To>
statslib_inline
ullint_t
spacing(const BlazeMat<eT,To>& X)
{
    return X.spacing();
}
#endif

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, int iTr, int iTc>
statslib_inline
ullint_t
spacing(const EigenMat<eT,iTr,iTc>& X)
{
    return n_cols(X);
}
#endif
