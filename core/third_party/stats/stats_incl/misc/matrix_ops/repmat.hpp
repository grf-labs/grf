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
// matrix repmat

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT>
statslib_inline
ArmaMat<eT>
repmat(const ArmaMat<eT>& X, ullint_t N, ullint_t K)
{
    return arma::repmat(X,N,K);
}
#endif

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, bool To>
statslib_inline
BlazeMat<eT,To>
repmat(const BlazeMat<eT,To>& X, ullint_t N, ullint_t K)
{
    ullint_t X_rows = X.rows();
    ullint_t X_cols = X.columns();

    BlazeMat<eT,To> mat_out(X_rows*N,X_cols*K);

    for (ullint_t j=ullint_t(0); j < N; ++j) {
        submatrix( mat_out, j*X_rows, 0U, X_rows,  X_cols ) = X;
    }

    if (K > ullint_t(1)) {
        for (ullint_t j=ullint_t(0); j < K; ++j) {
            submatrix( mat_out, 0U, j*X_cols, X_rows*N,  X_cols ) = submatrix( mat_out, 0U, X_cols, X_rows*N,  X_cols );
        }
    }

    return mat_out;
}
#endif

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, int iTr, int iTc>
statslib_inline
EigenMat<eT,iTr,iTc>
repmat(const EigenMat<eT,iTr,iTc>& X, ullint_t N, ullint_t K)
{
    return X.replicate(N,K);
}
#endif
