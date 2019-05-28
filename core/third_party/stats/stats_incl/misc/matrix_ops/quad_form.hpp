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
// x'Ax or x'(A^{-1})x, where x is (n x 1) and A is (n x n)

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT>
statslib_inline
double
quad_form(const ArmaMat<eT>& x, const ArmaMat<eT>& A, const bool inv_A = true)
{
    if (inv_A) {
        return static_cast<double>( arma::dot(x,arma::solve(A,x)) );
    } else {
        return static_cast<double>( arma::dot(x,A*x) );
    }
}
#endif

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename vT, typename mT>
statslib_inline
double
quad_form(const vT& x, const mT& A, const bool inv_A = true)
{
    if (inv_A) {
        // return static_cast<double>( blaze::dot(x, blaze::inv(A)*x) );
        mT res = blaze::trans(x) * blaze::inv(A) * x;
        return static_cast<double>( res(0,0) );
    } else {
        mT res = blaze::trans(x) * A * x;
        return static_cast<double>( res(0,0) );
    }
}
#endif

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename vT, typename mT>
statslib_inline
double
quad_form(const vT& x, const mT& A, const bool inv_A = true)
{
    if (inv_A) {
        const vT Ax = A.colPivHouseholderQr().solve(x);
        return static_cast<double>( (x.transpose() * Ax)(0) );
    } else {
        return static_cast<double>( (x.transpose() * A * x)(0) );
    }
}
#endif
