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
 * cdf of the Binomial distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
pbinom_compute(const llint_t x, const llint_t n_trials_par, const T prob_par, const llint_t count)
noexcept
{
    return( count == x ? \
                dbinom(count,n_trials_par,prob_par,false) : 
                dbinom(count,n_trials_par,prob_par,false) + pbinom_compute(x,n_trials_par,prob_par,count+1) );
}

template<typename T>
statslib_constexpr
T
pbinom_vals_check(const llint_t x, const llint_t n_trials_par, const T prob_par, const bool log_form)
noexcept
{
    return( !binom_sanity_check(n_trials_par,prob_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            x < llint_t(0) ? \
                log_zero_if<T>(log_form) :
            x >= n_trials_par ? \
                log_if(T(1),log_form) : // includes pbinom(0,0,.) case
            //
            n_trials_par == llint_t(1) ? \
                pbern(x,prob_par,log_form) :
            //
            log_if(pbinom_compute(x,n_trials_par,prob_par,llint_t(0)), log_form) );
}

}

/**
 * @brief Distribution function of the Binomial distribution
 *
 * @param x a real-valued input.
 * @param n_trials_par the number of trials, a non-negative integral-valued input.
 * @param prob_par the probability parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return the cumulative distribution function evaluated at \c x.
 * 
 * Example:
 * \code{.cpp} stats::pbinom(2,4,0.4,false); \endcode
 */

template<typename T>
statslib_constexpr
T
pbinom(const llint_t x, const llint_t n_trials_par, const T prob_par, const bool log_form)
noexcept
{
    return internal::pbinom_vals_check(x,n_trials_par,prob_par,log_form);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename rT>
statslib_inline
void
pbinom_vec(const eT* __stats_pointer_settings__ vals_in, const llint_t n_trials_par, const T1 prob_par, const bool log_form, 
                 rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(pbinom,vals_in,vals_out,num_elem,n_trials_par,prob_par,log_form);
}
#endif

}

/**
 * @brief Distribution function of the Binomial distribution
 *
 * @param x a standard vector.
 * @param n_trials_par the number of trials, a non-negative integral-valued input.
 * @param prob_par the probability parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return a vector of CDF values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<int> x = {2, 3, 4};
 * stats::pbinom(x,5,0.5,false);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
std::vector<rT>
pbinom(const std::vector<eT>& x, const llint_t n_trials_par, const T1 prob_par, const bool log_form)
{
    STDVEC_DIST_FN(pbinom_vec,n_trials_par,prob_par,log_form);
}
#endif

/**
 * @brief Distribution function of the Binomial distribution
 *
 * @param X a matrix of input values.
 * @param n_trials_par the number of trials, a non-negative integral-valued input.
 * @param prob_par the probability parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return a matrix of CDF values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::pbinom(X,5,0.5,false);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
ArmaMat<rT>
pbinom(const ArmaMat<eT>& X, const llint_t n_trials_par, const T1 prob_par, const bool log_form)
{
    ARMA_DIST_FN(pbinom_vec,n_trials_par,prob_par,log_form);
}

template<typename mT, typename tT, typename T1>
statslib_inline
mT
pbinom(const ArmaGen<mT,tT>& X, const llint_t n_trials_par, const T1 prob_par, const bool log_form)
{
    return pbinom(X.eval(),n_trials_par,prob_par,log_form);
}
#endif

/**
 * @brief Distribution function of the Binomial distribution
 *
 * @param X a matrix of input values.
 * @param n_trials_par the number of trials, a non-negative integral-valued input.
 * @param prob_par the probability parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return a matrix of CDF values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::pbinom(X,5,0.5,false);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
pbinom(const BlazeMat<eT,To>& X, const llint_t n_trials_par, const T1 prob_par, const bool log_form)
{
    BLAZE_DIST_FN(pbinom_vec,n_trials_par,prob_par,log_form);
}
#endif

/**
 * @brief Distribution function of the Binomial distribution
 *
 * @param X a matrix of input values.
 * @param n_trials_par the number of trials, a non-negative integral-valued input.
 * @param prob_par the probability parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return a matrix of CDF values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::pbinom(X,5,0.5,false);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
pbinom(const EigenMat<eT,iTr,iTc>& X, const llint_t n_trials_par, const T1 prob_par, const bool log_form)
{
    EIGEN_DIST_FN(pbinom_vec,n_trials_par,prob_par,log_form);
}
#endif
