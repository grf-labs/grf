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
 * Sample from an inverse-Wishart distribution
 */

/**
 * @brief Random sampling function for the Inverse-Wishart distribution
 *
 * @param Psi_par a positive semi-definite scale matrix.
 * @param nu_par the degrees of parameter, a real-valued input.
 * @param pre_chol indicate whether \c Psi_par is passed in lower triangular (Cholesky) format.
 *
 * @return a pseudo-random draw from the Inverse-Wishart distribution.
 */

template<typename mT, typename pT, typename not_arma_mat<mT>::type*>
statslib_inline
mT
rinvwish(const mT& Psi_par, const pT nu_par, const bool pre_chol)
{
    typedef return_t<pT> eT;
    const ullint_t K = mat_ops::n_rows(Psi_par);
    
    mT chol_Psi_inv;
    if (pre_chol) {
        chol_Psi_inv = Psi_par; // should be lower triangular
    } else {
        chol_Psi_inv = mat_ops::chol(mat_ops::inv(Psi_par)); // will be lower triangular
    }

    //

    rand_engine_t engine(std::random_device{}());

    mT A;
    mat_ops::zeros(A,K,K);

    for (ullint_t i=1U; i < K; i++) {
        for (ullint_t j=0U; j < i; j++) {
            A(i,j) = rnorm<eT>(eT(0),eT(1),engine);
        }
    }
    
    for (ullint_t i=0U; i < K; i++) {
        A(i,i) = std::sqrt(rchisq<eT>(eT(nu_par-i),engine));
    }

    chol_Psi_inv = chol_Psi_inv*A;

    //

    mT mat_out_inv = chol_Psi_inv * mat_ops::trans(chol_Psi_inv); // avoid Glue issues
    
    return mat_ops::inv( mat_out_inv );
}

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename mT, typename eT, typename pT>
statslib_inline
mT
rinvwish(const ArmaMat<eT>& Psi_par, const pT nu_par, const bool pre_chol)
{
    const ullint_t K = Psi_par.n_rows;
    
    ArmaMat<eT> chol_Psi = (pre_chol) ? Psi_par : arma::chol(arma::inv(Psi_par),"lower"); // should be lower-triangular

    //

    rand_engine_t engine(std::random_device{}());

    ArmaMat<eT> A = arma::zeros(K,K);

    for (ullint_t i=1U; i < K; i++) {
        for (ullint_t j=0U; j < i; j++) {
            A(i,j) = rnorm<eT>(eT(0),eT(1),engine);
        }
    }
    
    for (ullint_t i=0U; i < K; i++) {
        A(i,i) = std::sqrt(rchisq<eT>(eT(nu_par-i),engine));
    }

    chol_Psi = chol_Psi*A;

    //
    
    return arma::inv(chol_Psi * chol_Psi.t());
}
#endif
