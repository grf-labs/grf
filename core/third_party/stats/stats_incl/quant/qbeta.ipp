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
 * quantile function of the beta distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
qbeta_compute(const T p, const T a_par, const T b_par)
noexcept
{
    return gcem::incomplete_beta_inv(a_par,b_par,p);
}

template<typename T>
statslib_constexpr
T
qbeta_limit_vals(const T p, const T a_par, const T b_par)
noexcept
{
    return( a_par == T(0) && b_par == T(0) ? \
                p < T(0.5) ? \
                    T(0) : 
                p > T(0.5) ? \
                    T(1) :
                // p == T(0.5)
                    T(0.5) :
            //
            a_par == T(0) || (GCINT::is_posinf(b_par) && !GCINT::is_posinf(a_par)) ? \
                T(0) :
            //
            b_par == T(0) || (GCINT::is_posinf(a_par) && !GCINT::is_posinf(b_par)) ? \
                T(1) :
            // a, b == +Inf
                T(0.5) );
}

template<typename T>
statslib_constexpr
T
qbeta_vals_check(const T p, const T a_par, const T b_par)
noexcept
{
    return( !beta_sanity_check(a_par,b_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            !prob_val_check(p) ? \
                STLIM<T>::quiet_NaN() :
            //
            p == T(0) ? \
                T(0) :
            p == T(1) ? \
                T(1) :
            //
            (a_par == T(0) || b_par == T(0) || GCINT::any_posinf(a_par,b_par)) ? \
                qbeta_limit_vals(p,a_par,b_par) :
            //
            qbeta_compute(p,a_par,b_par) );
}

template<typename T1, typename T2, typename T3, typename TC = common_return_t<T1,T2,T3>>
statslib_constexpr
TC
qbeta_type_check(const T1 p, const T2 a_par, const T3 b_par)
noexcept
{
    return qbeta_vals_check(static_cast<TC>(p),static_cast<TC>(a_par),static_cast<TC>(b_par));
}

}

/**
 * @brief Quantile function of the Beta distribution
 *
 * @param p a real-valued input.
 * @param a_par shape parameter, a real-valued input.
 * @param b_par shape parameter, a real-valued input.
 *
 * @return the quantile function evaluated at \c p.
 * 
 * Example:
 * \code{.cpp} stats::qbeta(0.5,3.0,2.0); \endcode
 */

template<typename T1, typename T2, typename T3>
statslib_constexpr
common_return_t<T1,T2,T3>
qbeta(const T1 p, const T2 a_par, const T3 b_par)
noexcept
{
    return internal::qbeta_type_check(p,a_par,b_par);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
void
qbeta_vec(const eT* __stats_pointer_settings__ vals_in, const T1 a_par, const T2 b_par, 
                rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(qbeta,vals_in,vals_out,num_elem,a_par,b_par);
}
#endif

}

/**
 * @brief Quantile function of the Beta distribution
 *
 * @param x a standard vector.
 * @param a_par a real-valued shape parameter.
 * @param b_par a real-valued shape parameter.
 *
 * @return a vector of quantile values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<double> x = {0.3, 0.5, 0.9};
 * stats::qbeta(x,3.0,2.0);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
std::vector<rT>
qbeta(const std::vector<eT>& x, const T1 a_par, const T2 b_par)
{
    STDVEC_DIST_FN(qbeta_vec,a_par,b_par);
}
#endif

/**
 * @brief Quantile function of the Beta distribution
 *
 * @param X a matrix of input values.
 * @param a_par a real-valued shape parameter.
 * @param b_par a real-valued shape parameter.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * arma::mat X = { {0.2,  0.7,  0.1},
 *                 {0.9,  0.3,  0.87} };
 * stats::qbeta(X,3.0,2.0);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
ArmaMat<rT>
qbeta(const ArmaMat<eT>& X, const T1 a_par, const T2 b_par)
{
    ARMA_DIST_FN(qbeta_vec,a_par,b_par);
}

template<typename mT, typename tT, typename T1, typename T2>
statslib_inline
mT
qbeta(const ArmaGen<mT,tT>& X, const T1 a_par, const T2 b_par)
{
    return qbeta(X.eval(),a_par,b_par);
}
#endif

/**
 * @brief Quantile function of the Beta distribution
 *
 * @param X a matrix of input values.
 * @param a_par a real-valued shape parameter.
 * @param b_par a real-valued shape parameter.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::qbeta(X,3.0,2.0);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
qbeta(const BlazeMat<eT,To>& X, const T1 a_par, const T2 b_par)
{
    BLAZE_DIST_FN(qbeta_vec,a_par,b_par);
}
#endif

/**
 * @brief Quantile function of the Beta distribution
 *
 * @param X a matrix of input values.
 * @param a_par a real-valued shape parameter.
 * @param b_par a real-valued shape parameter.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::qbeta(X,3.0,2.0);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
qbeta(const EigenMat<eT,iTr,iTc>& X, const T1 a_par, const T2 b_par)
{
    EIGEN_DIST_FN(qbeta_vec,a_par,b_par);
}
#endif
