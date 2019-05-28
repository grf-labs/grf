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
 * quantile function of the gamma distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
qgamma_compute(const T p, const T shape_par, const T scale_par)
noexcept
{
    return( scale_par*gcem::incomplete_gamma_inv(shape_par,p) );
}

template<typename T>
statslib_constexpr
T
qgamma_limit_vals(const T shape_par, const T scale_par)
noexcept
{
    return( // here: 0 < p < 1
            shape_par == T(0) ? \
                GCINT::is_posinf(scale_par) ? \
                    STLIM<T>::quiet_NaN() :
                    T(0) :
            // shape or scale == +Inf
                STLIM<T>::quiet_NaN() );
}

template<typename T>
statslib_constexpr
T
qgamma_vals_check(const T p, const T shape_par, const T scale_par)
noexcept
{
    return( !gamma_sanity_check(shape_par,scale_par) ? \
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
            GCINT::any_posinf(shape_par,scale_par) || shape_par == T(0) ? \
                qgamma_limit_vals(shape_par,scale_par) :
            //
            qgamma_compute(p,shape_par,scale_par) );
}

template<typename T1, typename T2, typename T3, typename TC = common_return_t<T1,T2,T3>>
statslib_constexpr
TC
qgamma_type_check(const T1 x, const T2 shape_par, const T3 scale_par)
noexcept
{
    return qgamma_vals_check(static_cast<TC>(x),static_cast<TC>(shape_par),static_cast<TC>(scale_par));
}

}

/**
 * @brief Quantile function of the Gamma distribution
 *
 * @param p a real-valued input.
 * @param shape_par the shape parameter, a real-valued input.
 * @param scale_par the scale parameter, a real-valued input.
 *
 * @return the quantile function evaluated at \c p.
 *
 * Example:
 * \code{.cpp} stats::qgamma(0.4,2,3); \endcode
 */

template<typename T1, typename T2, typename T3>
statslib_constexpr
common_return_t<T1,T2,T3>
qgamma(const T1 p, const T2 shape_par, const T3 scale_par)
noexcept
{
    return internal::qgamma_type_check(p,shape_par,scale_par);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
void
qgamma_vec(const eT* __stats_pointer_settings__ vals_in, const T1 shape_par, const T2 scale_par, 
                 rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(qgamma,vals_in,vals_out,num_elem,shape_par,scale_par);
}
#endif

}

/**
 * @brief Quantile function of the Gamma distribution
 *
 * @param x a standard vector.
 * @param shape_par the shape parameter, a real-valued input.
 * @param scale_par the scale parameter, a real-valued input.
 *
 * @return a vector of quantile values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<double> x = {0.3, 0.5, 0.9};
 * stats::qgamma(x,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
std::vector<rT>
qgamma(const std::vector<eT>& x, const T1 shape_par, const T2 scale_par)
{
    STDVEC_DIST_FN(qgamma_vec,shape_par,scale_par);
}
#endif

/**
 * @brief Quantile function of the Gamma distribution
 *
 * @param X a matrix of input values.
 * @param shape_par the shape parameter, a real-valued input.
 * @param scale_par the scale parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * arma::mat X = { {0.2,  0.7,  0.1},
 *                 {0.9,  0.3,  0.87} };
 * stats::qgamma(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
ArmaMat<rT>
qgamma(const ArmaMat<eT>& X, const T1 shape_par, const T2 scale_par)
{
    ARMA_DIST_FN(qgamma_vec,shape_par,scale_par);
}

template<typename mT, typename tT, typename T1, typename T2>
statslib_inline
mT
qgamma(const ArmaGen<mT,tT>& X, const T1 shape_par, const T2 scale_par)
{
    return qgamma(X.eval(),shape_par,scale_par);
}
#endif

/**
 * @brief Quantile function of the Gamma distribution
 *
 * @param X a matrix of input values.
 * @param shape_par the shape parameter, a real-valued input.
 * @param scale_par the scale parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::qgamma(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
qgamma(const BlazeMat<eT,To>& X, const T1 shape_par, const T2 scale_par)
{
    BLAZE_DIST_FN(qgamma_vec,shape_par,scale_par);
}
#endif

/**
 * @brief Quantile function of the Gamma distribution
 *
 * @param X a matrix of input values.
 * @param shape_par the shape parameter, a real-valued input.
 * @param scale_par the scale parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::qgamma(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
qgamma(const EigenMat<eT,iTr,iTc>& X, const T1 shape_par, const T2 scale_par)
{
    EIGEN_DIST_FN(qgamma_vec,shape_par,scale_par);
}
#endif
