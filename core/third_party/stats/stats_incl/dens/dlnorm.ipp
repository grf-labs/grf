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
 * pdf of the univariate log-normal distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
dlnorm_log_check(const T x, const T mu_par, const T sigma_par, const bool log_form)
noexcept
{
    return( log_form == true ? \
                dnorm(stmath::log(x),mu_par,sigma_par,true) - stmath::log(x) :
                dnorm(stmath::log(x),mu_par,sigma_par,false) / x );
}

template<typename T>
statslib_constexpr
T
dlnorm_limit_vals(const T log_x, const T mu_par, const T sigma_par)
noexcept
{
    return( // sigma == Inf
            GCINT::is_posinf(sigma_par) ? \
                T(0) :
            // sigma == 0 and log(x) - mu == 0
            sigma_par == T(0) ? \
                GCINT::all_posinf(log_x,mu_par) || GCINT::all_neginf(log_x,mu_par) ? \
                    STLIM<T>::infinity() :
                    T(0) :
            // sigma finite and > 0; log(x) == mu == Inf
            GCINT::all_posinf(log_x,mu_par) ? \
                STLIM<T>::quiet_NaN() :
            //
                T(0) );
}

template<typename T>
statslib_constexpr
T
dlnorm_vals_check(const T x, const T mu_par, const T sigma_par, const bool log_form)
noexcept
{
    return( !lnorm_sanity_check(x,mu_par,sigma_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            x < T(0) ? \
                log_zero_if<T>(log_form) :
            //
            GCINT::any_inf(stmath::log(x),mu_par,sigma_par) || sigma_par == T(0) ? \
                log_if(dlnorm_limit_vals(stmath::log(x),mu_par,sigma_par),log_form) :
            //
            dlnorm_log_check(x,mu_par,sigma_par,log_form) );
}

template<typename T1, typename T2, typename T3, typename TC = common_return_t<T1,T2,T3>>
statslib_constexpr
TC
dlnorm_type_check(const T1 x, const T2 mu_par, const T3 sigma_par, const bool log_form)
noexcept
{
    return dlnorm_vals_check(static_cast<TC>(x),static_cast<TC>(mu_par),
                             static_cast<TC>(sigma_par),log_form);
}

}

/**
 * @brief Density function of the Log-Normal distribution
 *
 * @param x a real-valued input.
 * @param mu_par the mean parameter, a real-valued input.
 * @param sigma_par the standard deviation parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return the density function evaluated at \c x.
 * 
 * Example:
 * \code{.cpp} stats::dlnorm(2.0,1.0,2.0,false); \endcode
 */

template<typename T1, typename T2, typename T3>
statslib_constexpr
common_return_t<T1,T2,T3>
dlnorm(const T1 x, const T2 mu_par, const T3 sigma_par, const bool log_form)
noexcept
{
    return internal::dlnorm_type_check(x,mu_par,sigma_par,log_form);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
void
dlnorm_vec(const eT* __stats_pointer_settings__ vals_in, const T1 mu_par, const T2 sigma_par, const bool log_form, 
                 rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(dlnorm,vals_in,vals_out,num_elem,mu_par,sigma_par,log_form);
}
#endif

}

/**
 * @brief Density function of the Log-Normal distribution
 *
 * @param x a standard vector.
 * @param mu_par the mean parameter, a real-valued input.
 * @param sigma_par the standard deviation parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a vector of density function values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<double> x = {0.0, 1.0, 2.0};
 * stats::dlnorm(x,1.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
std::vector<rT>
dlnorm(const std::vector<eT>& x, const T1 mu_par, const T2 sigma_par, const bool log_form)
{
    STDVEC_DIST_FN(dlnorm_vec,mu_par,sigma_par,log_form);
}
#endif

/**
 * @brief Density function of the Log-Normal distribution
 *
 * @param X a matrix of input values.
 * @param mu_par the mean parameter, a real-valued input.
 * @param sigma_par the standard deviation parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * arma::mat X = { {0.2, 1.7, 0.1},
 *                 {0.9, 4.0, 0.3} };
 * stats::dlnorm(X,1.0,1.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
ArmaMat<rT>
dlnorm(const ArmaMat<eT>& X, const T1 mu_par, const T2 sigma_par, const bool log_form)
{
    ARMA_DIST_FN(dlnorm_vec,mu_par,sigma_par,log_form);
}

template<typename mT, typename tT, typename T1, typename T2>
statslib_inline
mT
dlnorm(const ArmaGen<mT,tT>& X, const T1 mu_par, const T2 sigma_par, const bool log_form)
{
    return dlnorm(X.eval(),mu_par,sigma_par,log_form);
}
#endif

/**
 * @brief Density function of the Log-Normal distribution
 *
 * @param X a matrix of input values.
 * @param mu_par the mean parameter, a real-valued input.
 * @param sigma_par the standard deviation parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::dlnorm(X,1.0,1.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
dlnorm(const BlazeMat<eT,To>& X, const T1 mu_par, const T2 sigma_par, const bool log_form)
{
    BLAZE_DIST_FN(dlnorm,mu_par,sigma_par,log_form);
}
#endif

/**
 * @brief Density function of the Log-Normal distribution
 *
 * @param X a matrix of input values.
 * @param mu_par the mean parameter, a real-valued input.
 * @param sigma_par the standard deviation parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::dlnorm(X,1.0,1.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
dlnorm(const EigenMat<eT,iTr,iTc>& X, const T1 mu_par, const T2 sigma_par, const bool log_form)
{
    EIGEN_DIST_FN(dlnorm_vec,mu_par,sigma_par,log_form);
}
#endif
