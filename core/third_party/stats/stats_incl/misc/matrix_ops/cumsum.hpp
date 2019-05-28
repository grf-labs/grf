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
// Cumulative sum (assumes a vector form)

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT>
statslib_inline
std::vector<eT>
cumsum(const std::vector<eT>& X)
{
    std::vector<eT> mat_out = X;

    eT* mem_out = mat_out.data();

    for (ullint_t j=ullint_t(1); j < X.rows()*X.columns(); ++j)
    {
        mem_out[j] += mem_out[j-1];
    }

    return mat_out;
}
#endif

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT>
statslib_inline
ArmaMat<eT>
cumsum(const ArmaMat<eT>& X)
{
    return arma::cumsum(X);
}
#endif

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, bool To>
statslib_inline
BlazeMat<eT,To>
cumsum(const BlazeMat<eT,To>& X)
{
    BlazeMat<eT,To> mat_out = X;

    eT* mem_out = mat_out.data();

    for (ullint_t j=ullint_t(1); j < X.rows()*X.columns(); ++j)
    {
        mem_out[j] += mem_out[j-1];
    }

    return mat_out;
}
#endif

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, int iTr, int iTc>
statslib_inline
EigenMat<eT,iTr,iTc>
cumsum(const EigenMat<eT,iTr,iTc>& X)
{
    EigenMat<eT,iTr,iTc> mat_out = X;

    eT* mem_out = mat_out.data();

    for (ullint_t j=ullint_t(1); j < X.rows()*X.cols(); ++j)
    {
        mem_out[j] += mem_out[j-1];
    }

    return mat_out;
}
#endif
