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
 * cdf of the Weibull distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
pweibull_compute(const T x, const T shape_par)
noexcept
{
    return( - stmath::expm1(-stmath::pow(x,shape_par)) );
}

template<typename T>
statslib_constexpr
T
pweibull_limit_vals(const T x, const T shape_par, const T scale_par)
noexcept
{
    return( GCINT::is_posinf(x) ? \
                GCINT::is_posinf(scale_par) ? \
                    STLIM<T>::quiet_NaN() :
                    T(1) :
            // x finite but > 0
            GCINT::is_posinf(scale_par) ? \
                T(0) :
            //
            GCINT::is_posinf(shape_par) ? \
                x / scale_par > T(1) ? \
                    T(1) :
                x / scale_par == T(1) ? \
                    - stmath::expm1(-1.0) :
                    T(0) :
            //
            T(0) );
}

template<typename T>
statslib_constexpr
T
pweibull_vals_check(const T x, const T shape_par, const T scale_par, const bool log_form)
noexcept
{
    return( !weibull_sanity_check(x,shape_par,scale_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            x < STLIM<T>::epsilon() ? \
                log_zero_if<T>(log_form) :
            //
            GCINT::any_inf(x,shape_par,scale_par) ? \
                log_if(pweibull_limit_vals(x,shape_par,scale_par), log_form) :
            //
            log_if(pweibull_compute(x/scale_par,shape_par), log_form) );
}

template<typename T1, typename T2, typename T3, typename TC = common_return_t<T1,T2,T3>>
statslib_constexpr
TC
pweibull_type_check(const T1 x, const T2 shape_par, const T3 scale_par, const bool log_form)
noexcept
{
    return pweibull_vals_check(static_cast<TC>(x),static_cast<TC>(shape_par),
                               static_cast<TC>(scale_par),log_form);
}

}

/**
 * @brief Distribution function of the Weibull distribution
 *
 * @param x a real-valued input.
 * @param shape_par the shape parameter, a real-valued input.
 * @param scale_par the scale parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return the cumulative distribution function evaluated at \c x.
 * 
 * Example:
 * \code{.cpp} stats::pweibull(1.0,2.0,3.0,false); \endcode
 */

template<typename T1, typename T2, typename T3>
statslib_constexpr
common_return_t<T1,T2,T3>
pweibull(const T1 x, const T2 shape_par, const T3 scale_par, const bool log_form)
noexcept
{
    return internal::pweibull_type_check(x,shape_par,scale_par,log_form);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
void
pweibull_vec(const eT* __stats_pointer_settings__ vals_in, const T1 shape_par, const T2 scale_par, const bool log_form, 
                   rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(pweibull,vals_in,vals_out,num_elem,shape_par,scale_par,log_form);
}
#endif

}

/**
 * @brief Distribution function of the Weibull distribution
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
 * stats::pweibull(x,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
std::vector<rT>
pweibull(const std::vector<eT>& x, const T1 shape_par, const T2 scale_par, const bool log_form)
{
    STDVEC_DIST_FN(pweibull_vec,shape_par,scale_par,log_form);
}
#endif

/**
 * @brief Distribution function of the Weibull distribution
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
 * stats::pweibull(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
ArmaMat<rT>
pweibull(const ArmaMat<eT>& X, const T1 shape_par, const T2 scale_par, const bool log_form)
{
    ARMA_DIST_FN(pweibull_vec,shape_par,scale_par,log_form);
}

template<typename mT, typename tT, typename T1, typename T2>
statslib_inline
mT
pweibull(const ArmaGen<mT,tT>& X, const T1 shape_par, const T2 scale_par, const bool log_form)
{
    return pweibull(X.eval(),shape_par,scale_par,log_form);
}
#endif

/**
 * @brief Distribution function of the Weibull distribution
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
 * stats::pweibull(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
pweibull(const BlazeMat<eT,To>& X, const T1 shape_par, const T2 scale_par, const bool log_form)
{
    BLAZE_DIST_FN(pweibull_vec,shape_par,scale_par,log_form);
}
#endif

/**
 * @brief Distribution function of the Weibull distribution
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
 * stats::pweibull(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
pweibull(const EigenMat<eT,iTr,iTc>& X, const T1 shape_par, const T2 scale_par, const bool log_form)
{
    EIGEN_DIST_FN(pweibull_vec,shape_par,scale_par,log_form);
}
#endif
