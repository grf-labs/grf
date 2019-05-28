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
// log-determinant

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT>
statslib_inline
eT
log_det(const ArmaMat<eT>& X)
{
    ArmaMat<eT> vec = mat_ops::chol(X).diag();
    return 2.0*arma::accu(arma::log(vec));
}
#endif

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, bool To>
statslib_inline
eT
log_det(const BlazeMat<eT,To>& X)
{
    BlazeMat<eT,To> chol_sig = mat_ops::chol(X);
    auto vec = blaze::diagonal(chol_sig);
    eT total = eT(0);
    for (auto it=vec.cbegin(); it!=vec.cend(); ++it)
    {
        total += std::log(*it);
    }
    return 2.0*total;
}
#endif

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, int iTr, int iTc>
statslib_inline
eT
log_det(const EigenMat<eT,iTr,iTc>& X)
{
    return ( mat_ops::chol(X).diagonal().array().log()*2.0 ).sum();
}
#endif
