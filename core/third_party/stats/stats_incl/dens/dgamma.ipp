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
 * pdf of the gamma distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
dgamma_log_compute(const T x, const T shape_par, const T scale_par)
noexcept
{
    return( - stmath::lgamma(shape_par) - shape_par*stmath::log(scale_par) \
            + (shape_par-T(1))*stmath::log(x) - x/scale_par );
}

template<typename T>
statslib_constexpr
T
dgamma_limit_vals(const T x, const T shape_par, const T scale_par)
noexcept
{
    return( shape_par == T(0) ? \
                x == T(0) ? \
                    STLIM<T>::infinity() :
                    T(0) :
            // x == 0
            shape_par < T(1) ? \
                STLIM<T>::infinity() :
            shape_par == T(1) ? \
                T(1) / scale_par :
            // x == 0 and shape_par > 1
                T(0) );
}

template<typename T>
statslib_constexpr
T
dgamma_vals_check(const T x, const T shape_par, const T scale_par, const bool log_form)
noexcept
{
    return( !gamma_sanity_check(x,shape_par,scale_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            x < T(0) ? \
                log_zero_if<T>(log_form) :
            //
            x == T(0) || shape_par == T(0) ? \
                log_if(dgamma_limit_vals(x,shape_par,scale_par), log_form) :
            //
            GCINT::any_posinf(x,shape_par,scale_par) ? \
                log_zero_if<T>(log_form) :
            //
            exp_if(dgamma_log_compute(x,shape_par,scale_par), !log_form) );
}

template<typename T1, typename T2, typename T3, typename TC = common_return_t<T1,T2,T3>>
statslib_constexpr
TC
dgamma_type_check(const T1 x, const T2 shape_par, const T3 scale_par, const bool log_form)
noexcept
{
    return dgamma_vals_check(static_cast<TC>(x),static_cast<TC>(shape_par),
                             static_cast<TC>(scale_par),log_form);
}

}

/**
 * @brief Density function of the Gamma distribution
 *
 * @param x a real-valued input.
 * @param shape_par the shape parameter, a real-valued input.
 * @param scale_par the scale parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return the density function evaluated at \c x.
 *
 * Example:
 * \code{.cpp} stats::dgamma(2,2,3,false); \endcode
 */

template<typename T1, typename T2, typename T3>
statslib_constexpr
common_return_t<T1,T2,T3>
dgamma(const T1 x, const T2 shape_par, const T3 scale_par, const bool log_form)
noexcept
{
    return internal::dgamma_type_check(x,shape_par,scale_par,log_form);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
void
dgamma_vec(const eT* __stats_pointer_settings__ vals_in, const T1 shape_par, const T2 scale_par, const bool log_form, 
                 rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(dgamma,vals_in,vals_out,num_elem,shape_par,scale_par,log_form);
}
#endif

}

/**
 * @brief Density function of the Gamma distribution
 *
 * @param x a standard vector.
 * @param shape_par the shape parameter, a real-valued input.
 * @param scale_par the scale parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a vector of density function values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<double> x = {1.8, 0.7, 4.2};
 * stats::dgamma(x,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
std::vector<rT>
dgamma(const std::vector<eT>& x, const T1 shape_par, const T2 scale_par, const bool log_form)
{
    STDVEC_DIST_FN(dgamma_vec,shape_par,scale_par,log_form);
}
#endif

/**
 * @brief Density function of the Gamma distribution
 *
 * @param X a matrix of input values.
 * @param shape_par the shape parameter, a real-valued input.
 * @param scale_par the scale parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * arma::mat X = { {1.8, 0.7, 4.2},
 *                 {0.3, 5.3, 3.7} };
 * stats::dgamma(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
ArmaMat<rT>
dgamma(const ArmaMat<eT>& X, const T1 shape_par, const T2 scale_par, const bool log_form)
{
    ARMA_DIST_FN(dgamma_vec,shape_par,scale_par,log_form);
}

template<typename mT, typename tT, typename T1, typename T2>
statslib_inline
mT
dgamma(const ArmaGen<mT,tT>& X, const T1 shape_par, const T2 scale_par, const bool log_form)
{
    return dgamma(X.eval(),shape_par,scale_par,log_form);
}
#endif

/**
 * @brief Density function of the Gamma distribution
 *
 * @param X a matrix of input values.
 * @param shape_par the shape parameter, a real-valued input.
 * @param scale_par the scale parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::dgamma(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
dgamma(const BlazeMat<eT,To>& X, const T1 shape_par, const T2 scale_par, const bool log_form)
{
    BLAZE_DIST_FN(dgamma,shape_par,scale_par,log_form);
}
#endif

/**
 * @brief Density function of the Gamma distribution
 *
 * @param X a matrix of input values.
 * @param shape_par the shape parameter, a real-valued input.
 * @param scale_par the scale parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::dgamma(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
dgamma(const EigenMat<eT,iTr,iTc>& X, const T1 shape_par, const T2 scale_par, const bool log_form)
{
    EIGEN_DIST_FN(dgamma_vec,shape_par,scale_par,log_form);
}
#endif
