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
 * pdf of the Laplace distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
dlaplace_log_compute(const T x, const T mu_par, const T sigma_par)
noexcept
{
    return( - stmath::log(2*sigma_par) - stmath::abs(x - mu_par) / sigma_par );
}

template<typename T>
statslib_constexpr
T
dlaplace_limit_vals(const T x, const T mu_par, const T sigma_par)
noexcept
{
    return( // sigma == Inf
            GCINT::is_posinf(sigma_par) ? \
                T(0) :
            // sigma finite; x == mu == Inf or -Inf 
            GCINT::all_posinf(x,mu_par) || GCINT::all_neginf(x,mu_par) ? \
                STLIM<T>::quiet_NaN() :
            // sigma finite; x and mu have opposite Inf signs
            // GCINT::is_posinf(x) || GCINT::is_neginf(mu_par) ?
                // T(0) :
            // GCINT::is_neginf(x) || GCINT::is_posinf(mu_par) ?
                // T(0) :
            // sigma == 0 and x-mu == 0
            sigma_par == T(0) && x == mu_par ? \
                STLIM<T>::infinity() :
            //
                T(0) );
}

template<typename T>
statslib_constexpr
T
dlaplace_vals_check(const T x, const T mu_par, const T sigma_par, const bool log_form)
noexcept
{
    return( !laplace_sanity_check(x,mu_par,sigma_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            GCINT::any_inf(x,mu_par,sigma_par) || sigma_par == T(0) ? \
                log_if(dlaplace_limit_vals(x,mu_par,sigma_par),log_form) :
            //
            exp_if(dlaplace_log_compute(x,mu_par,sigma_par), !log_form) );
}

template<typename T1, typename T2, typename T3, typename TC = common_return_t<T1,T2,T3>>
statslib_constexpr
TC
dlaplace_type_check(const T1 x, const T2 mu_par, const T3 sigma_par, const bool log_form)
noexcept
{
    return dlaplace_vals_check(static_cast<TC>(x),static_cast<TC>(mu_par),
                               static_cast<TC>(sigma_par),log_form);
}

}

/**
 * @brief Density function of the Laplace distribution
 *
 * @param x a real-valued input.
 * @param mu_par the location parameter, a real-valued input.
 * @param sigma_par the scale parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return the density function evaluated at \c x.
 * 
 * Example:
 * \code{.cpp} stats::dlaplace(0.7,1.0,2.0,false); \endcode
 */

template<typename T1, typename T2, typename T3>
statslib_constexpr
common_return_t<T1,T2,T3>
dlaplace(const T1 x, const T2 mu_par, const T3 sigma_par, const bool log_form)
noexcept
{
    return internal::dlaplace_type_check(x,mu_par,sigma_par,log_form);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
void
dlaplace_vec(const eT* __stats_pointer_settings__ vals_in, const T1 mu_par, const T2 sigma_par, const bool log_form, 
                   rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(dlaplace,vals_in,vals_out,num_elem,mu_par,sigma_par,log_form);
}
#endif

}

/**
 * @brief Density function of the Laplace distribution
 *
 * @param x a standard vector.
 * @param mu_par the location parameter, a real-valued input.
 * @param sigma_par the scale parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a vector of density function values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<double> x = {0.0, 1.0, 2.0};
 * stats::dlaplace(x,1.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
std::vector<rT>
dlaplace(const std::vector<eT>& x, const T1 mu_par, const T2 sigma_par, const bool log_form)
{
    STDVEC_DIST_FN(dlaplace_vec,mu_par,sigma_par,log_form);
}
#endif

/**
 * @brief Density function of the Laplace distribution
 *
 * @param X a matrix of input values.
 * @param mu_par the location parameter, a real-valued input.
 * @param sigma_par the scale parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * arma::mat X = { {0.2, -1.7,  0.1},
 *                 {0.9,  4.0, -0.3} };
 * stats::dlaplace(X,1.0,1.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
ArmaMat<rT>
dlaplace(const ArmaMat<eT>& X, const T1 mu_par, const T2 sigma_par, const bool log_form)
{
    ARMA_DIST_FN(dlaplace_vec,mu_par,sigma_par,log_form);
}

template<typename mT, typename tT, typename T1, typename T2>
statslib_inline
mT
dlaplace(const ArmaGen<mT,tT>& X, const T1 mu_par, const T2 sigma_par, const bool log_form)
{
    return dlaplace(X.eval(),mu_par,sigma_par,log_form);
}
#endif

/**
 * @brief Density function of the Laplace distribution
 *
 * @param X a matrix of input values.
 * @param mu_par the location parameter, a real-valued input.
 * @param sigma_par the scale parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::dlaplace(X,1.0,1.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
dlaplace(const BlazeMat<eT,To>& X, const T1 mu_par, const T2 sigma_par, const bool log_form)
{
    BLAZE_DIST_FN(dlaplace,mu_par,sigma_par,log_form);
}
#endif

/**
 * @brief Density function of the Laplace distribution
 *
 * @param X a matrix of input values.
 * @param mu_par the location parameter, a real-valued input.
 * @param sigma_par the scale parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::dlaplace(X,1.0,1.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
dlaplace(const EigenMat<eT,iTr,iTc>& X, const T1 mu_par, const T2 sigma_par, const bool log_form)
{
    EIGEN_DIST_FN(dlaplace_vec,mu_par,sigma_par,log_form);
}
#endif
