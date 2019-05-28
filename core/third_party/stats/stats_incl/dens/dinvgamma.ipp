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
 * pdf of the inverse-gamma distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
dinvgamma_log_compute(const T x, const T shape_par, const T rate_par)
noexcept
{
    return( - stmath::lgamma(shape_par) + shape_par*stmath::log(rate_par) \
            + (-shape_par-T(1))*stmath::log(x) - rate_par/x );
}

template<typename T>
statslib_constexpr
T
dinvgamma_limit_vals(const T x, const T shape_par, const T rate_par)
noexcept
{
    return( shape_par == T(0) ? \
                x > T(0) ? \
                    STLIM<T>::infinity() :
                    T(0) :
            // x = 0
            rate_par == T(0) ? \
                STLIM<T>::infinity() :
                T(0) );
}

template<typename T>
statslib_constexpr
T
dinvgamma_vals_check(const T x, const T shape_par, const T rate_par, const bool log_form)
noexcept
{
    return( !invgamma_sanity_check(x,shape_par,rate_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            x < T(0) ? \
                log_zero_if<T>(log_form) :
            //
            x == T(0) || shape_par == T(0) ? \
                log_if(dinvgamma_limit_vals(x,shape_par,rate_par), log_form) :
            //
            rate_par == T(0) ? \
                log_zero_if<T>(log_form) :
            //
            GCINT::any_posinf(x,shape_par,rate_par) ? \
                log_zero_if<T>(log_form) :
            //
            exp_if(dinvgamma_log_compute(x,shape_par,rate_par), !log_form) );
}

template<typename T1, typename T2, typename T3, typename TC = common_return_t<T1,T2,T3>>
statslib_constexpr
TC
dinvgamma_type_check(const T1 x, const T2 shape_par, const T3 rate_par, const bool log_form)
noexcept
{
    return dinvgamma_vals_check(static_cast<TC>(x),static_cast<TC>(shape_par),
                               static_cast<TC>(rate_par),log_form);
}

}

/**
 * @brief Density function of the Inverse-Gamma distribution
 *
 * @param x a real-valued input.
 * @param shape_par the shape parameter, a real-valued input.
 * @param rate_par the rate parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return the density function evaluated at \c x.
 *
 * Example:
 * \code{.cpp} stats::dinvgamma(1.5,2,1,false); \endcode
 */

template<typename T1, typename T2, typename T3>
statslib_constexpr
common_return_t<T1,T2,T3>
dinvgamma(const T1 x, const T2 shape_par, const T3 rate_par, const bool log_form)
noexcept
{
    return internal::dinvgamma_type_check(x,shape_par,rate_par,log_form);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
void
dinvgamma_vec(const eT* __stats_pointer_settings__ vals_in, const T1 shape_par, const T2 rate_par, const bool log_form, 
                    rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(dinvgamma,vals_in,vals_out,num_elem,shape_par,rate_par,log_form);
}
#endif

}

/**
 * @brief Density function of the Inverse-Gamma distribution
 *
 * @param x a standard vector.
 * @param shape_par the shape parameter, a real-valued input.
 * @param rate_par the rate parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a vector of density function values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<double> x = {1.8, 0.7, 4.2};
 * stats::dinvgamma(x,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
std::vector<rT>
dinvgamma(const std::vector<eT>& x, const T1 shape_par, const T2 rate_par, const bool log_form)
{
    STDVEC_DIST_FN(dinvgamma_vec,shape_par,rate_par,log_form);
}
#endif

/**
 * @brief Density function of the Inverse-Gamma distribution
 *
 * @param X a matrix of input values.
 * @param shape_par the shape parameter, a real-valued input.
 * @param rate_par the rate parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * arma::mat X = { {1.8, 0.7, 4.2},
 *                 {0.3, 5.3, 3.7} };
 * stats::dinvgamma(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
ArmaMat<rT>
dinvgamma(const ArmaMat<eT>& X, const T1 shape_par, const T2 rate_par, const bool log_form)
{
    ARMA_DIST_FN(dinvgamma_vec,shape_par,rate_par,log_form);
}

template<typename mT, typename tT, typename T1, typename T2>
statslib_inline
mT
dinvgamma(const ArmaGen<mT,tT>& X, const T1 shape_par, const T2 rate_par, const bool log_form)
{
    return dinvgamma(X.eval(),shape_par,rate_par,log_form);
}
#endif

/**
 * @brief Density function of the Inverse-Gamma distribution
 *
 * @param X a matrix of input values.
 * @param shape_par the shape parameter, a real-valued input.
 * @param rate_par the rate parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::dinvgamma(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
dinvgamma(const BlazeMat<eT,To>& X, const T1 shape_par, const T2 rate_par, const bool log_form)
{
    BLAZE_DIST_FN(dinvgamma,shape_par,rate_par,log_form);
}
#endif

/**
 * @brief Density function of the Inverse-Gamma distribution
 *
 * @param X a matrix of input values.
 * @param shape_par the shape parameter, a real-valued input.
 * @param rate_par the rate parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::dinvgamma(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
dinvgamma(const EigenMat<eT,iTr,iTc>& X, const T1 shape_par, const T2 rate_par, const bool log_form)
{
    EIGEN_DIST_FN(dinvgamma_vec,shape_par,rate_par,log_form);
}
#endif
