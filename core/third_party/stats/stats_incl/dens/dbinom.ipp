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
 * pdf of the Binomial distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
dbinom_log_compute(const llint_t x, const llint_t n_trials_par, const T prob_par)
noexcept
{
    return( x == llint_t(0) ? \
                    n_trials_par * stmath::log(T(1.0) - prob_par) :
            //
            x == n_trials_par ? \
                x * stmath::log(prob_par) :
            //
            gcem::log_binomial_coef(n_trials_par,x) + x*stmath::log(prob_par) \
                + (n_trials_par - x)*stmath::log(T(1) - prob_par) );
}

template<typename T>
statslib_constexpr
T
dbinom_limit_vals(const llint_t x)
noexcept
{
    return( x == llint_t(0) ? \
                T(1) :
                T(0) );
}

template<typename T>
statslib_constexpr
T
dbinom_vals_check(const llint_t x, const llint_t n_trials_par, const T prob_par, const bool log_form)
noexcept
{
    return( !binom_sanity_check(n_trials_par,prob_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            GCINT::is_nan(x) ? \
                STLIM<T>::quiet_NaN() :
            //
            x < llint_t(0) || x > n_trials_par ? \
                log_zero_if<T>(log_form) :
            //
            n_trials_par == llint_t(0) ? \
                log_if(dbinom_limit_vals<T>(x),log_form) :
            //
            n_trials_par == llint_t(1) ? \
                dbern(x,prob_par,log_form) :
            //
            exp_if(dbinom_log_compute(x,n_trials_par,prob_par), !log_form) );
}

}

/**
 * @brief Density function of the Binomial distribution
 *
 * @param x a real-valued input.
 * @param n_trials_par the number of trials, a non-negative integral-valued input.
 * @param prob_par the probability parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return the density function evaluated at \c x.
 * 
 * Example:
 * \code{.cpp} stats::dbinom(2,4,0.4,false); \endcode
 */

template<typename T>
statslib_constexpr
return_t<T>
dbinom(const llint_t x, const llint_t n_trials_par, const T prob_par, const bool log_form)
noexcept
{
    return internal::dbinom_vals_check(x,n_trials_par,static_cast<return_t<T>>(prob_par),log_form);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename rT>
statslib_inline
void
dbinom_vec(const eT* __stats_pointer_settings__ vals_in, const llint_t n_trials_par, const T1 prob_par, const bool log_form, 
                 rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(dbinom,vals_in,vals_out,num_elem,n_trials_par,prob_par,log_form);
}
#endif

}

/**
 * @brief Density function of the Binomial distribution
 *
 * @param x a standard vector.
 * @param n_trials_par the number of trials, a non-negative integral-valued input.
 * @param prob_par the probability parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a vector of density function values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<int> x = {2, 3, 4};
 * stats::dbinom(x,5,0.5,false);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
std::vector<rT>
dbinom(const std::vector<eT>& x, const llint_t n_trials_par, const T1 prob_par, const bool log_form)
{
    STDVEC_DIST_FN(dbinom_vec,n_trials_par,prob_par,log_form);
}
#endif

/**
 * @brief Density function of the Binomial distribution
 *
 * @param X a matrix of input values.
 * @param n_trials_par the number of trials, a non-negative integral-valued input.
 * @param prob_par the probability parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::dbinom(X,5,0.5,false);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
ArmaMat<rT>
dbinom(const ArmaMat<eT>& X, const llint_t n_trials_par, const T1 prob_par, const bool log_form)
{
    ARMA_DIST_FN(dbinom_vec,n_trials_par,prob_par,log_form);
}

template<typename mT, typename tT, typename T1>
statslib_inline
mT
dbinom(const ArmaGen<mT,tT>& X, const llint_t n_trials_par, const T1 prob_par, const bool log_form)
{
    return dbinom(X.eval(),n_trials_par,prob_par,log_form);
}
#endif

/**
 * @brief Density function of the Binomial distribution
 *
 * @param X a matrix of input values.
 * @param n_trials_par the number of trials, a non-negative integral-valued input.
 * @param prob_par the probability parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::dbinom(X,5,0.5,false);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
dbinom(const BlazeMat<eT,To>& X, const llint_t n_trials_par, const T1 prob_par, const bool log_form)
{
    BLAZE_DIST_FN(dbinom,n_trials_par,prob_par,log_form);
}
#endif

/**
 * @brief Density function of the Binomial distribution
 *
 * @param X a matrix of input values.
 * @param n_trials_par the number of trials, a non-negative integral-valued input.
 * @param prob_par the probability parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::dbinom(X,5,0.5,false);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
dbinom(const EigenMat<eT,iTr,iTc>& X, const llint_t n_trials_par, const T1 prob_par, const bool log_form)
{
    EIGEN_DIST_FN(dbinom_vec,n_trials_par,prob_par,log_form);
}
#endif
