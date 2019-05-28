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
 * pdf of the chi-squared distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
dchisq_compute(const T x, const T dof_par)
noexcept
{
    return( - stmath::lgamma(0.5*dof_par) - T(0.5)*dof_par*T(GCEM_LOG_2) \
                + (T(0.5)*dof_par - T(1))*stmath::log(x) - x / T(2.0) );
}

template<typename T>
statslib_constexpr
T
dchisq_limit_vals(const T dof_par)
noexcept
{
    return( // x == 0 cases 
            dof_par < T(2) ? \
                STLIM<T>::infinity() :
            dof_par == T(2) ? \
                T(0.5) :
                T(0) );
}

template<typename T>
statslib_constexpr
T
dchisq_vals_check(const T x, const T dof_par, const bool log_form)
noexcept
{
    return( !chisq_sanity_check(x,dof_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            x < T(0) ? \
                log_zero_if<T>(log_form) :
            //
            x == T(0) ? \
                log_if(dchisq_limit_vals(dof_par), log_form) :
            // dof == +Inf or x == +Inf
            GCINT::any_posinf(x,dof_par) ? \
                log_zero_if<T>(log_form) :
            //
            exp_if(dchisq_compute(x,dof_par), !log_form) );
}

template<typename T1, typename T2, typename TC = common_return_t<T1,T2>>
statslib_constexpr
TC
dchisq_type_check(const T1 x, const T2 dof_par, const bool log_form)
noexcept
{
    return dchisq_vals_check(static_cast<TC>(x),static_cast<TC>(dof_par),log_form);
}

}

/**
 * @brief Density function of the Chi-squared distribution
 *
 * @param x a real-valued input.
 * @param dof_par the degrees of freedom parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return the density function evaluated at \c x.
 * 
 * Example:
 * \code{.cpp} stats::dchisq(4,5,false); \endcode
 */

template<typename T1, typename T2>
statslib_constexpr
common_return_t<T1,T2>
dchisq(const T1 x, const T2 dof_par, const bool log_form)
noexcept
{
    return internal::dchisq_type_check(x,dof_par,log_form);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename rT>
statslib_inline
void
dchisq_vec(const eT* __stats_pointer_settings__ vals_in, const T1 dof_par, const bool log_form, 
                 rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(dchisq,vals_in,vals_out,num_elem,dof_par,log_form);
}
#endif

}

/**
 * @brief Density function of the Chi-squared distribution
 *
 * @param x a standard vector.
 * @param dof_par the degrees of freedom parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a vector of density function values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<double> x = {1.8, 0.7, 4.2};
 * stats::dchisq(x,4,false);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
std::vector<rT>
dchisq(const std::vector<eT>& x, const T1 dof_par, const bool log_form)
{
    STDVEC_DIST_FN(dchisq_vec,dof_par,log_form);
}
#endif

/**
 * @brief Density function of the Chi-squared distribution
 *
 * @param X a matrix of input values.
 * @param dof_par the degrees of freedom parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * arma::mat X = { {1.8, 0.7, 4.2},
 *                 {0.3, 5.3, 3.7} };
 * stats::dchisq(X,4,false);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
ArmaMat<rT>
dchisq(const ArmaMat<eT>& X, const T1 dof_par, const bool log_form)
{
    ARMA_DIST_FN(dchisq_vec,dof_par,log_form);
}

template<typename mT, typename tT, typename T1>
statslib_inline
mT
dchisq(const ArmaGen<mT,tT>& X, const T1 dof_par, const bool log_form)
{
    return dchisq(X.eval(),dof_par,log_form);
}
#endif

/**
 * @brief Density function of the Chi-squared distribution
 *
 * @param X a matrix of input values.
 * @param dof_par the degrees of freedom parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::dchisq(X,4,false);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
dchisq(const BlazeMat<eT,To>& X, const T1 dof_par, const bool log_form)
{
    BLAZE_DIST_FN(dchisq,dof_par,log_form);
}
#endif

/**
 * @brief Density function of the Chi-squared distribution
 *
 * @param X a matrix of input values.
 * @param dof_par the degrees of freedom parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::dchisq(X,4,false);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
dchisq(const EigenMat<eT,iTr,iTc>& X, const T1 dof_par, const bool log_form)
{
    EIGEN_DIST_FN(dchisq_vec,dof_par,log_form);
}
#endif
