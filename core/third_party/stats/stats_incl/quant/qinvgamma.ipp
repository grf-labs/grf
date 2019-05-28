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
 * quantile function of the inverse-gamma distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
qinvgamma_compute(const T p, const T shape_par, const T rate_par)
noexcept
{
    return( rate_par / gcem::incomplete_gamma_inv(shape_par,T(1)-p) );
}

template<typename T>
statslib_constexpr
T
qinvgamma_limit_vals(const T shape_par, const T rate_par)
noexcept
{
    return( // here: 0 < p < 1
            rate_par == T(0) ? \
                shape_par == T(0) ? \
                    STLIM<T>::quiet_NaN() :
                    T(0) :
            // shape == 0; shape or rate == +Inf
                STLIM<T>::quiet_NaN() );
}

template<typename T>
statslib_constexpr
T
qinvgamma_vals_check(const T p, const T shape_par, const T rate_par)
noexcept
{
    return( !invgamma_sanity_check(shape_par,rate_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            !prob_val_check(p) ? \
                STLIM<T>::quiet_NaN() :
            //
            p == T(0) ? \
                T(0) :
            p == T(1) ? \
                STLIM<T>::infinity() :
            //
            GCINT::any_posinf(shape_par,rate_par) || shape_par == T(0) || rate_par == T(0) ? \
                qinvgamma_limit_vals(shape_par,rate_par) :
            //
            qinvgamma_compute(p,shape_par,rate_par) );
}

template<typename T1, typename T2, typename T3, typename TC = common_return_t<T1,T2,T3>>
statslib_constexpr
TC
qinvgamma_type_check(const T1 x, const T2 shape_par, const T3 rate_par)
noexcept
{
    return qinvgamma_vals_check(static_cast<TC>(x),static_cast<TC>(shape_par),static_cast<TC>(rate_par));
}

}

/**
 * @brief Quantile function of the Inverse-Gamma distribution
 *
 * @param p a real-valued input.
 * @param shape_par the shape parameter, a real-valued input.
 * @param rate_par the rate parameter, a real-valued input.
 *
 * @return the quantile function evaluated at \c p.
 *
 * Example:
 * \code{.cpp} stats::qinvgamma(0.5,2,1); \endcode
 */

template<typename T1, typename T2, typename T3>
statslib_constexpr
common_return_t<T1,T2,T3>
qinvgamma(const T1 p, const T2 shape_par, const T3 rate_par)
noexcept
{
    return internal::qinvgamma_type_check(p,shape_par,rate_par);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
void
qinvgamma_vec(const eT* __stats_pointer_settings__ vals_in, const T1 shape_par, const T2 rate_par, 
                    rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(qinvgamma,vals_in,vals_out,num_elem,shape_par,rate_par);
}
#endif

}

/**
 * @brief Quantile function of the Inverse-Gamma distribution
 *
 * @param x a standard vector.
 * @param shape_par the shape parameter, a real-valued input.
 * @param rate_par the rate parameter, a real-valued input.
 *
 * @return a vector of quantile values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<double> x = {0.3, 0.5, 0.9};
 * stats::qinvgamma(x,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
std::vector<rT>
qinvgamma(const std::vector<eT>& x, const T1 shape_par, const T2 rate_par)
{
    STDVEC_DIST_FN(qinvgamma_vec,shape_par,rate_par);
}
#endif

/**
 * @brief Quantile function of the Inverse-Gamma distribution
 *
 * @param X a matrix of input values.
 * @param shape_par the shape parameter, a real-valued input.
 * @param rate_par the rate parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * arma::mat X = { {0.2,  0.7,  0.1},
 *                 {0.9,  0.3,  0.87} };
 * stats::qinvgamma(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
ArmaMat<rT>
qinvgamma(const ArmaMat<eT>& X, const T1 shape_par, const T2 rate_par)
{
    ARMA_DIST_FN(qinvgamma_vec,shape_par,rate_par);
}

template<typename mT, typename tT, typename T1, typename T2>
statslib_inline
mT
qinvgamma(const ArmaGen<mT,tT>& X, const T1 shape_par, const T2 rate_par)
{
    return qinvgamma(X.eval(),shape_par,rate_par);
}
#endif

/**
 * @brief Quantile function of the Inverse-Gamma distribution
 *
 * @param X a matrix of input values.
 * @param shape_par the shape parameter, a real-valued input.
 * @param rate_par the rate parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::qinvgamma(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
qinvgamma(const BlazeMat<eT,To>& X, const T1 shape_par, const T2 rate_par)
{
    BLAZE_DIST_FN(qinvgamma_vec,shape_par,rate_par);
}
#endif

/**
 * @brief Quantile function of the Inverse-Gamma distribution
 *
 * @param X a matrix of input values.
 * @param shape_par the shape parameter, a real-valued input.
 * @param rate_par the rate parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::qinvgamma(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
qinvgamma(const EigenMat<eT,iTr,iTc>& X, const T1 shape_par, const T2 rate_par)
{
    EIGEN_DIST_FN(qinvgamma_vec,shape_par,rate_par);
}
#endif
