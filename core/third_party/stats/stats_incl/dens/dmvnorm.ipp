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
 * pdf of the Multivariate Normal distribution
 */

/**
 * @brief Density function of the Multivariate-Normal distribution
 *
 * @param X a column vector.
 * @param mu_par mean vector.
 * @param Sigma_par the covariance matrix.
 * @param log_form return the log-density or the true form.
 *
 * @return the density function evaluated at \c X.
 */

template<typename vT, typename mT, typename eT>
statslib_inline
eT
dmvnorm(const vT& X, const vT& mu_par, const mT& Sigma_par, bool log_form)
{
    const ullint_t K = mat_ops::n_rows(X);

    //

    const eT cons_term = static_cast<eT>( -eT(0.5)*K*GCEM_LOG_2PI );
    const vT X_cent = X - mu_par; // avoids issues like Mat vs eGlue in templates
    const eT quad_term = mat_ops::quad_form(X_cent, Sigma_par, true);
    
    eT ret = cons_term - eT(0.5) * ( mat_ops::log_det(Sigma_par) + quad_term );

    if (!log_form) {
        ret = std::exp(ret);
        
        if (std::isinf(ret)) {
            ret = std::numeric_limits<eT>::max();
        }
    }

    //
    
    return ret;
}
