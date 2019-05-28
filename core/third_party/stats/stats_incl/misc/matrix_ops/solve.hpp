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
// matrix 'solve': A*X = B

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT>
statslib_inline
ArmaMat<eT>
solve(const ArmaMat<eT>& A, const ArmaMat<eT>& B)
{
    return arma::solve(A,B);
}
#endif

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, bool To>
statslib_inline
BlazeMat<eT,To>
solve(const BlazeMat<eT,To>& A, const BlazeMat<eT,To>& B)
{
    return blaze::inv(A)*B;
}
#endif

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
// template<typename eT, int iTr, int iTc>
// statslib_inline
// EigenMat<eT,iTr,iTc>
// solve(const EigenMat<eT,iTr,iTr>& A, const EigenMat<eT,iTr,iTc>& B)
// {
//     // note A must be square, so iTr=iTc for A
//     return A.colPivHouseholderQr().solve(B);
// }

template<typename mT, typename vT>
statslib_inline
vT
solve(const mT& A, const vT& B)
{
    // note A must be square, so iTr=iTc for A
    return A.colPivHouseholderQr().solve(B);
}
#endif
