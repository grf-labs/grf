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
// element-wise exponential

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT>
statslib_inline
std::vector<eT>
exp(const std::vector<eT>& X)
{
    std::vector<eT> vec_out = X;
    std::for_each(vec_out.begin(), vec_out.end(), [](eT& x){ x = std::exp(x);});
    return vec_out;
}
#endif

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT>
statslib_inline
ArmaMat<eT>
exp(const ArmaMat<eT>& X)
{
    return arma::exp(X);
}
#endif

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, bool To>
statslib_inline
BlazeMat<eT,To>
exp(const BlazeMat<eT,To>& X)
{
    return blaze::exp(X);
}
#endif

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, int iTr, int iTc>
statslib_inline
EigenMat<eT,iTr,iTc>
exp(const EigenMat<eT,iTr,iTc>& X)
{
    return X.array().exp().matrix();
}
#endif
