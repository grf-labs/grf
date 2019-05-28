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
 * cdf of the gamma distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
pgamma_compute(const T x, const T shape_par, const T scale_par)
noexcept
{
    return gcem::incomplete_gamma(shape_par,x/scale_par);
}

template<typename T>
statslib_constexpr
T
pgamma_limit_vals(const T x, const T shape_par, const T scale_par)
noexcept
{
    return( // here: x > 0
            shape_par == T(0) ? \
                GCINT::is_posinf(scale_par) ? \
                    T(0) :
                    T(1) :
            // x == +Inf, pars finite
            GCINT::is_posinf(x) && GCINT::all_finite(shape_par,scale_par) ? \
                T(1) :
            // x, shape_par == +Inf, scale finite
            GCINT::all_posinf(x,shape_par) ? \
                T(1) :
            // shape or scale == +Inf; x finite
                T(0) );
}

template<typename T>
statslib_constexpr
T
pgamma_vals_check(const T x, const T shape_par, const T scale_par, const bool log_form)
noexcept
{
    return( !gamma_sanity_check(x,shape_par,scale_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            GCINT::all_posinf(x,scale_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            STLIM<T>::epsilon() > x ? \
                log_zero_if<T>(log_form) :
            //
            GCINT::any_posinf(x,shape_par,scale_par) || shape_par == T(0) ? \
                log_if(pgamma_limit_vals(x,shape_par,scale_par), log_form) :
            //
            log_if(pgamma_compute(x,shape_par,scale_par), log_form) );
}

template<typename T1, typename T2, typename T3, typename TC = common_return_t<T1,T2,T3>>
statslib_constexpr
TC
pgamma_type_check(const T1 x, const T2 shape_par, const T3 scale_par, const bool log_form)
noexcept
{
    return pgamma_vals_check(static_cast<TC>(x),static_cast<TC>(shape_par),
                             static_cast<TC>(scale_par),log_form);
}

}

/**
 * @brief Distribution function of the Gamma distribution
 *
 * @param x a real-valued input.
 * @param shape_par the shape parameter, a real-valued input.
 * @param scale_par the scale parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return the cumulative distribution function evaluated at \c x.
 *
 * Example:
 * \code{.cpp} stats::pgamma(2,2,3,false); \endcode
 */

template<typename T1, typename T2, typename T3>
statslib_constexpr
common_return_t<T1,T2,T3>
pgamma(const T1 x, const T2 shape_par, const T3 scale_par, const bool log_form)
noexcept
{
    return internal::pgamma_type_check(x,shape_par,scale_par,log_form);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
void
pgamma_vec(const eT* __stats_pointer_settings__ vals_in, const T1 shape_par, const T2 scale_par, const bool log_form, 
                 rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(pgamma,vals_in,vals_out,num_elem,shape_par,scale_par,log_form);
}
#endif

}

/**
 * @brief Distribution function of the Gamma distribution
 *
 * @param x a standard vector.
 * @param shape_par the shape parameter, a real-valued input.
 * @param scale_par the scale parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return a vector of CDF values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<double> x = {1.8, 0.7, 4.2};
 * stats::pgamma(x,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
std::vector<rT>
pgamma(const std::vector<eT>& x, const T1 shape_par, const T2 scale_par, const bool log_form)
{
    STDVEC_DIST_FN(pgamma_vec,shape_par,scale_par,log_form);
}
#endif

/**
 * @brief Distribution function of the Gamma distribution
 *
 * @param X a matrix of input values.
 * @param shape_par the shape parameter, a real-valued input.
 * @param scale_par the scale parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return a matrix of CDF values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * arma::mat X = { {1.8, 0.7, 4.2},
 *                 {0.3, 5.3, 3.7} };
 * stats::pgamma(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
ArmaMat<rT>
pgamma(const ArmaMat<eT>& X, const T1 shape_par, const T2 scale_par, const bool log_form)
{
    ARMA_DIST_FN(pgamma_vec,shape_par,scale_par,log_form);
}

template<typename mT, typename tT, typename T1, typename T2>
statslib_inline
mT
pgamma(const ArmaGen<mT,tT>& X, const T1 shape_par, const T2 scale_par, const bool log_form)
{
    return pgamma(X.eval(),shape_par,scale_par,log_form);
}
#endif

/**
 * @brief Distribution function of the Gamma distribution
 *
 * @param X a matrix of input values.
 * @param shape_par the shape parameter, a real-valued input.
 * @param scale_par the scale parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return a matrix of CDF values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::pgamma(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
pgamma(const BlazeMat<eT,To>& X, const T1 shape_par, const T2 scale_par, const bool log_form)
{
    BLAZE_DIST_FN(pgamma_vec,shape_par,scale_par,log_form);
}
#endif

/**
 * @brief Distribution function of the Gamma distribution
 *
 * @param X a matrix of input values.
 * @param shape_par the shape parameter, a real-valued input.
 * @param scale_par the scale parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return a matrix of CDF values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::pgamma(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
pgamma(const EigenMat<eT,iTr,iTc>& X, const T1 shape_par, const T2 scale_par, const bool log_form)
{
    EIGEN_DIST_FN(pgamma_vec,shape_par,scale_par,log_form);
}
#endif
