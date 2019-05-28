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
// vector variance

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT>
statslib_inline
eT
var(const std::vector<eT>& X)
{
    eT mean_val = mean(X);
    eT sq_val = sqaccu(X) / static_cast<eT>(X.size());
    return sq_val - mean_val*mean_val;
}
#endif

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT>
statslib_inline
eT
var(const ArmaMat<eT>& X)
{
    return arma::as_scalar(arma::var(X));
}
#endif

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, bool To>
statslib_inline
eT
var(const BlazeMat<eT,To>& X)
{
    eT mean_val = mean(X);
    eT sq_val = sqaccu(X) / static_cast<eT>(n_elem(X));
    return sq_val - mean_val*mean_val;
}
#endif

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, int iTr, int iTc>
statslib_inline
eT
var(const EigenMat<eT,iTr,iTc>& X)
{
    eT mean_val = mean(X);
    eT sq_val = sqaccu(X) / static_cast<eT>(n_elem(X));
    return sq_val - mean_val*mean_val;
}
#endif
