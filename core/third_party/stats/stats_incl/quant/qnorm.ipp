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
 * quantile function of the univariate normal distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
qnorm_compute(const T p, const T mu_par, const T sigma_par)
noexcept
{
    return( mu_par + T(GCEM_SQRT_2)*sigma_par*gcem::erf_inv( T(2)*p - T(1) ) );
}

template<typename T>
statslib_constexpr
T
qnorm_vals_check(const T p, const T mu_par, const T sigma_par)
noexcept
{
    return( !logis_sanity_check(mu_par,sigma_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            !prob_val_check(p) ? \
                STLIM<T>::quiet_NaN() :
            //
            GCINT::is_posinf(sigma_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            sigma_par == T(0) ? \
                mu_par :
            //
            p == T(0) ? \
                GCINT::is_posinf(mu_par) ? \
                      STLIM<T>::quiet_NaN() :
                    - STLIM<T>::infinity() :
            p == T(1) ? \
                GCINT::is_neginf(mu_par) ? \
                    STLIM<T>::quiet_NaN() :
                    STLIM<T>::infinity() :
            //
            qnorm_compute(p,mu_par,sigma_par) );
}

template<typename T1, typename T2, typename T3, typename TC = common_return_t<T1,T2,T3>>
statslib_constexpr
TC
qnorm_type_check(const T1 p, const T2 mu_par, const T3 sigma_par)
noexcept
{
    return qnorm_vals_check(static_cast<TC>(p),static_cast<TC>(mu_par),static_cast<TC>(sigma_par));
}

}

/**
 * Quantile function of the Normal distribution
 *
 * @param p a real-valued input.
 * @param mu_par the mean parameter, a real-valued input.
 * @param sigma_par the standard deviation parameter, a real-valued input.
 *
 * @return the quantile function evaluated at \c p.
 *
 * Example:
 * \code{.cpp} stats::qnorm(0.5,1.0,2.0); \endcode
 */

template<typename T1, typename T2, typename T3>
statslib_constexpr
common_return_t<T1,T2,T3>
qnorm(const T1 p, const T2 mu_par, const T3 sigma_par)
noexcept
{
    return internal::qnorm_type_check(p,mu_par,sigma_par);
}

/**
 * Quantile function of the Standard Normal distribution
 *
 * @param p a real-valued input.
 *
 * @return the quantile function evaluated at \c p.
 *
 * Example:
 * \code{.cpp} stats::qnorm(0.5); \endcode
 */

template<typename T>
statslib_constexpr
return_t<T>
qnorm(const T p)
noexcept
{
    return qnorm(p,return_t<T>(0),return_t<T>(1));
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
void
qnorm_vec(const eT* __stats_pointer_settings__ vals_in, const T1 mu_par, const T2 sigma_par, 
                rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(qnorm,vals_in,vals_out,num_elem,mu_par,sigma_par);
}
#endif

}

/**
 * @brief Quantile function of the Normal distribution
 *
 * @param x a standard vector.
 * @param mu_par the mean parameter, a real-valued input.
 * @param sigma_par the standard deviation parameter, a real-valued input.
 *
 * @return a vector of quantile values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<double> x = {0.1, 0.3, 0.7};
 * stats::qnorm(x,1.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
std::vector<rT>
qnorm(const std::vector<eT>& x, const T1 mu_par, const T2 sigma_par)
{
    STDVEC_DIST_FN(qnorm_vec,mu_par,sigma_par);
}
#endif

/**
 * @brief Quantile function of the Normal distribution
 *
 * @param X a matrix of input values.
 * @param mu_par the mean parameter, a real-valued input.
 * @param sigma_par the standard deviation parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * arma::mat X = { {0.2, 0.7, 0.9},
 *                 {0.1, 0.8, 0.3} };
 * stats::qnorm(X,1.0,1.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
ArmaMat<rT>
qnorm(const ArmaMat<eT>& X, const T1 mu_par, const T2 sigma_par)
{
    ARMA_DIST_FN(qnorm_vec,mu_par,sigma_par);
}

template<typename mT, typename tT, typename T1, typename T2>
statslib_inline
mT
qnorm(const ArmaGen<mT,tT>& X, const T1 mu_par, const T2 sigma_par)
{
    return qnorm(X.eval(),mu_par,sigma_par);
}
#endif

/**
 * @brief Quantile function of the Normal distribution
 *
 * @param X a matrix of input values.
 * @param mu_par the mean parameter, a real-valued input.
 * @param sigma_par the standard deviation parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::qnorm(X,1.0,1.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
qnorm(const BlazeMat<eT,To>& X, const T1 mu_par, const T2 sigma_par)
{
    BLAZE_DIST_FN(qnorm_vec,mu_par,sigma_par);
}
#endif

/**
 * @brief Quantile function of the Normal distribution
 *
 * @param X a matrix of input values.
 * @param mu_par the mean parameter, a real-valued input.
 * @param sigma_par the standard deviation parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::qnorm(X,1.0,1.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
qnorm(const EigenMat<eT,iTr,iTc>& X, const T1 mu_par, const T2 sigma_par)
{
    EIGEN_DIST_FN(qnorm_vec,mu_par,sigma_par);
}
#endif
