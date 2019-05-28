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
 * cdf of the exponential distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
qexp_compute(const T p, const T rate_par)
noexcept
{
    return( - stmath::log( T(1) - p ) / rate_par );
}

template<typename T>
statslib_constexpr
T
qexp_vals_check(const T p, const T rate_par)
noexcept
{
    return( !exp_sanity_check(rate_par) ? \
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
            qexp_compute(p,rate_par) );
}

template<typename T1, typename T2, typename TC = common_return_t<T1,T2>>
statslib_constexpr
TC
qexp_type_check(const T1 p, const T2 rate_par)
noexcept
{
    return qexp_vals_check(static_cast<TC>(p),static_cast<TC>(rate_par));
}

}

/**
 * @brief Quantile function of the Exponential distribution
 *
 * @param p a real-valued input.
 * @param rate_par the rate parameter, a real-valued input.
 *
 * @return the quantile function evaluated at \c p.
 * 
 * Example:
 * \code{.cpp} stats::qexp(0.5,4.0); \endcode
 */

template<typename T1, typename T2>
statslib_constexpr
common_return_t<T1,T2>
qexp(const T1 p, const T2 rate_par)
noexcept
{
    return internal::qexp_type_check(p,rate_par);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename rT>
statslib_inline
void
qexp_vec(const eT* __stats_pointer_settings__ vals_in, const T1 rate_par, 
                 rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(qexp,vals_in,vals_out,num_elem,rate_par);
}
#endif

}

/**
 * @brief Quantile function of the Exponential distribution
 *
 * @param x a standard vector.
 * @param rate_par the rate parameter, a real-valued input.
 *
 * @return a vector of quantile values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<double> x = {0.3, 0.5, 0.8};
 * stats::qexp(x,4);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
std::vector<rT>
qexp(const std::vector<eT>& x, const T1 rate_par)
{
    STDVEC_DIST_FN(qexp_vec,rate_par);
}
#endif

/**
 * @brief Quantile function of the Exponential distribution
 *
 * @param X a matrix of input values.
 * @param rate_par the rate parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * arma::mat X = { {0.2, 0.7, 0.9},
 *                 {0.1, 0.8, 0.3} };
 * stats::qexp(X,4);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
ArmaMat<rT>
qexp(const ArmaMat<eT>& X, const T1 rate_par)
{
    ARMA_DIST_FN(qexp_vec,rate_par);
}

template<typename mT, typename tT, typename T1>
statslib_inline
mT
qexp(const ArmaGen<mT,tT>& X, const T1 rate_par)
{
    return qexp(X.eval(),rate_par);
}
#endif

/**
 * @brief Quantile function of the Exponential distribution
 *
 * @param X a matrix of input values.
 * @param rate_par the rate parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::qexp(X,4);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
qexp(const BlazeMat<eT,To>& X, const T1 rate_par)
{
    BLAZE_DIST_FN(qexp_vec,rate_par);
}
#endif

/**
 * @brief Quantile function of the Exponential distribution
 *
 * @param X a matrix of input values.
 * @param rate_par the rate parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::qexp(X,4);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
qexp(const EigenMat<eT,iTr,iTc>& X, const T1 rate_par)
{
    EIGEN_DIST_FN(qexp_vec,rate_par);
}
#endif
