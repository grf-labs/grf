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
// sum all elements and sum of squared values

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT>
statslib_inline
eT
accu(const std::vector<eT>& X)
{
    // const eT sum_val = std::accumulate(X.begin(), X.end(), eT(0));
    eT sum_val = eT(0);
    for (auto x : X)
        sum_val += x;
    return sum_val;
}

template<typename eT>
statslib_inline
eT
sqaccu(const std::vector<eT>& X)
{
    eT sum_val = eT(0);
    for (auto& x : X)
        sum_val += x*x;
    return sum_val;
}
#endif

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT>
statslib_inline
eT
accu(const ArmaMat<eT>& X)
{
    return arma::accu(X);
}

template<typename eT>
statslib_inline
eT
sqaccu(const ArmaMat<eT>& X)
{
    return arma::accu(arma::pow(X,2));
}
#endif

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, bool To>
statslib_inline
eT
accu(const BlazeMat<eT,To>& X)
{
    eT out_val = blaze::sum(X);
    return out_val;
}

template<typename eT, bool To>
statslib_inline
eT
sqaccu(const BlazeMat<eT,To>& X)
{
    eT out_val = blaze::sum(blaze::pow(X,2));
    return out_val;
}
#endif

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, int iTr, int iTc>
statslib_inline
eT
accu(const EigenMat<eT,iTr,iTc>& X)
{
    return X.sum();
}

template<typename eT, int iTr, int iTc>
statslib_inline
eT
sqaccu(const EigenMat<eT,iTr,iTc>& X)
{
    // const eT* vals = X.data();
    // eT out_val = eT(0);
    // for (ullint_t j=0U; j < n_elem(X); ++j) {
    //     out_val += vals[j]*vals[j];
    // }
    // return out_val;
    return (X.pow(2)).sum();
}
#endif
