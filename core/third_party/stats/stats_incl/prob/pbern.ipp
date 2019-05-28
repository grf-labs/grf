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
 * cdf of the Bernoulli distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
pbern_compute(const llint_t x, const T prob_par)
noexcept
{
    return( x < llint_t(0) ? \
                T(0) :
            //
            x >= llint_t(1) ? \
                T(1) : 
                T(1) - prob_par );
}

template<typename T>
statslib_constexpr
T
pbern_vals_check(const llint_t x, const T prob_par, const bool log_form)
noexcept
{
    return( !bern_sanity_check(prob_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            log_if(pbern_compute(x,prob_par), log_form) );
}

}

/**
 * @brief Distribution function of the Bernoulli distribution
 *
 * @param x a value equal to 0 or 1.
 * @param prob_par the probability parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return the cumulative distribution function evaluated at \c x.
 * 
 * Example:
 * \code{.cpp} stats::pbern(1,0.6,false); \endcode
 */

template<typename T>
statslib_constexpr
return_t<T>
pbern(const llint_t x, const T prob_par, const bool log_form)
noexcept
{
    return internal::pbern_vals_check(x,static_cast<return_t<T>>(prob_par),log_form);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename rT>
statslib_inline
void
pbern_vec(const eT* __stats_pointer_settings__ vals_in, const T1 prob_par, const bool log_form,
                rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(pbern,vals_in,vals_out,num_elem,prob_par,log_form);
}
#endif

}

/**
 * @brief Density function of the Bernoulli distribution
 *
 * @param x a standard vector.
 * @param prob_par the probability parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return a vector of CDF values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<int> x = {0, 1, 0};
 * stats::pbern(x,0.5,false);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
std::vector<rT>
pbern(const std::vector<eT>& x, const T1 prob_par, const bool log_form)
{
    STDVEC_DIST_FN(pbern_vec,prob_par,log_form);
}
#endif

/**
 * @brief Density function of the Bernoulli distribution
 *
 * @param X a matrix of input values.
 * @param prob_par the probability parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return a matrix of CDF values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * arma::mat X = { {1, 0, 1},
 *                 {0, 1, 0} };
 * stats::pbern(X,0.5,false);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
ArmaMat<rT>
pbern(const ArmaMat<eT>& X, const T1 prob_par, const bool log_form)
{
    ARMA_DIST_FN(pbern_vec,prob_par,log_form);
}

template<typename mT, typename tT, typename T1>
statslib_inline
mT
pbern(const ArmaGen<mT,tT>& X, const T1 prob_par, const bool log_form)
{
    return pbern(X.eval(),prob_par,log_form);
}
#endif

/**
 * @brief Density function of the Bernoulli distribution
 *
 * @param X a matrix of input values.
 * @param prob_par the probability parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return a matrix of CDF values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::pbern(X,0.5,false);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
pbern(const BlazeMat<eT,To>& X, const T1 prob_par, const bool log_form)
{
    BLAZE_DIST_FN(pbern_vec,prob_par,log_form);
}
#endif

/**
 * @brief Density function of the Bernoulli distribution
 *
 * @param X a matrix of input values.
 * @param prob_par the probability parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return a matrix of CDF values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::pbern(X,0.5,false);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
pbern(const EigenMat<eT,iTr,iTc>& X, const T1 prob_par, const bool log_form)
{
    EIGEN_DIST_FN(pbern_vec,prob_par,log_form);
}
#endif
