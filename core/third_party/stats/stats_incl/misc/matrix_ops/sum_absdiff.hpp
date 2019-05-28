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
// sum the absolute element-wise differences between two objects of the same dimensions

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT>
statslib_inline
eT
sum_absdiff(const std::vector<eT>& X, const std::vector<eT>& Y)
{
    eT val_out = eT(0);
    ullint_t n_elem = X.size(); // assumes dim(X) = dim(Y)

    for (ullint_t i=ullint_t(0); i < n_elem; ++i)
    {
        val_out += std::abs(X[i] - Y[i]);
    }

    return val_out;
}
#endif

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT>
statslib_inline
eT
sum_absdiff(const ArmaMat<eT>& X, const ArmaMat<eT>& Y)
{
    return arma::accu(arma::abs(X - Y));
}
#endif

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, bool To>
statslib_inline
eT
sum_absdiff(const BlazeMat<eT,To>& X, const BlazeMat<eT,To>& Y)
{
    return blaze::sum(blaze::abs(X-Y));
}
#endif

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, int iTr, int iTc>
statslib_inline
eT
sum_absdiff(const EigenMat<eT,iTr,iTc>& X, const EigenMat<eT,iTr,iTc>& Y)
{
    return (X - Y).array().abs().sum();
}
#endif
