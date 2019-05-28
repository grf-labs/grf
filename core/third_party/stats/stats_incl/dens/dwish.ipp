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
 * pdf of the Wishart distribution
 */

/**
 * @brief Density function of the Wishart distribution
 *
 * @param X a positive semi-definite matrix.
 * @param Psi_par a positive semi-definite scale matrix.
 * @param nu_par the degrees of parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return the density function evaluated at \c X.
 */

template<typename mT, typename pT, typename not_arma_mat<mT>::type*>
statslib_inline
return_t<pT>
dwish(const mT& X, const mT& Psi_par, const pT nu_par, const bool log_form)
{
    typedef return_t<pT> eT;

    const ullint_t K = mat_ops::n_rows(X);
    const eT nu_par_d2 = static_cast<eT>(nu_par) / eT(2);

    //

    const eT lmg_term = gcem::lmgamma(nu_par_d2, K);
    const eT norm_term = - nu_par_d2*std::log(mat_ops::det(Psi_par)) - nu_par_d2*K*GCEM_LOG_2 - lmg_term;

    eT ret = norm_term + eT(0.5) * ( (nu_par-K-1) * std::log(mat_ops::det(X)) - mat_ops::trace(mat_ops::solve(Psi_par,X)) );

    if (!log_form) {
        ret = std::exp(ret);
    }

    //

    return ret;
}

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename pT>
statslib_inline
eT
dwish(const ArmaMat<eT>& X, const ArmaMat<eT>& Psi_par, const pT nu_par, const bool log_form)
{
    const ullint_t K = X.n_rows;
    const eT nu_par_d2 = static_cast<eT>(nu_par) / eT(2);

    //

    const eT lmg_term = gcem::lmgamma(nu_par_d2, K);
    const eT norm_term = - nu_par_d2*std::log(arma::det(Psi_par)) - nu_par_d2*K*GCEM_LOG_2 - lmg_term;

    eT ret = norm_term + eT(0.5) * ( (nu_par-K-1) * std::log(arma::det(X)) - arma::trace(arma::solve(Psi_par,X)) );

    if (!log_form) {
        ret = std::exp(ret);
    }

    //

    return ret;
}
#endif
