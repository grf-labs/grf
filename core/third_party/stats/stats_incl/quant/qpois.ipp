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
 * quantile function of the Poisson distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
qpois_compute_right_search(const T p, const T rate_par, const T value, const llint_t count)
noexcept
{
    return( value <= p ? \
                qpois_compute_right_search(p,rate_par, ppois(count,rate_par,false), count + llint_t(1)) : 
            // else
            count > llint_t(0) ? 
                static_cast<T>(count - llint_t(1)) : 
                T(0) );
}

template<typename T>
statslib_constexpr
T
qpois_compute_search_begin(const T p, const T rate_par, const llint_t count)
noexcept
{
    return qpois_compute_right_search(p,rate_par, ppois(count,rate_par,false), count);
}

template<typename T>
statslib_constexpr
llint_t
qpois_trunc_normal_approx(const T p, const T rate_par)
noexcept
{
    return static_cast<llint_t>( stmath::max( T(0), qnorm(p,rate_par,stmath::sqrt(rate_par)) - T(3) ) );
}

template<typename T>
statslib_constexpr
T
qpois_vals_check(const T p, const T rate_par)
noexcept
{
    return( !pois_sanity_check(rate_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            !prob_val_check(p) ? \
                STLIM<T>::quiet_NaN() :
            //
            p == T(0) ? \
                T(0) :
            STLIM<T>::epsilon() > rate_par ? \
                T(0) :
            p == T(1) ? \
                STLIM<T>::infinity() :
            // rate < 11
            rate_par < T(11) ? \
                qpois_compute_right_search(p,rate_par,T(0),llint_t(0)) :
            // else use a normal approximation
                qpois_compute_search_begin(p,rate_par,qpois_trunc_normal_approx(p,rate_par)) );
}

template<typename T1, typename T2, typename TC = common_return_t<T1,T2>>
statslib_constexpr
TC
qpois_type_check(const T1 p, const T2 rate_par)
noexcept
{
    return qpois_vals_check(static_cast<TC>(p),static_cast<TC>(rate_par));
}

}

/**
 * @brief Quantile function of the Poisson distribution
 *
 * @param p a real-valued input.
 * @param rate_par the rate parameter, a real-valued input.
 *
 * @return the quantile function evaluated at \c p.
 * 
 * Example:
 * \code{.cpp} stats::qpois(0.6,10.0); \endcode
 */

template<typename T1, typename T2>
statslib_constexpr
common_return_t<T1,T2>
qpois(const T1 p, const T2 rate_par)
noexcept
{
    return internal::qpois_type_check(p,rate_par);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename rT>
statslib_inline
void
qpois_vec(const eT* __stats_pointer_settings__ vals_in, const T1 rate_par, 
                 rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(qpois,vals_in,vals_out,num_elem,rate_par);
}
#endif

}

/**
 * @brief Quantile function of the Poisson distribution
 *
 * @param x a standard vector.
 * @param rate_par the rate parameter, a real-valued input.
 *
 * @return a vector of quantile values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<double> x = {0.3, 0.5, 0.8};
 * stats::qpois(x,4);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
std::vector<rT>
qpois(const std::vector<eT>& x, const T1 rate_par)
{
    STDVEC_DIST_FN(qpois_vec,rate_par);
}
#endif

/**
 * @brief Quantile function of the Poisson distribution
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
 * stats::qpois(X,4);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
ArmaMat<rT>
qpois(const ArmaMat<eT>& X, const T1 rate_par)
{
    ARMA_DIST_FN(qpois_vec,rate_par);
}

template<typename mT, typename tT, typename T1>
statslib_inline
mT
qpois(const ArmaGen<mT,tT>& X, const T1 rate_par)
{
    return qpois(X.eval(),rate_par);
}
#endif

/**
 * @brief Quantile function of the Poisson distribution
 *
 * @param X a matrix of input values.
 * @param rate_par the rate parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::qpois(X,4);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
qpois(const BlazeMat<eT,To>& X, const T1 rate_par)
{
    BLAZE_DIST_FN(qpois_vec,rate_par);
}
#endif

/**
 * @brief Quantile function of the Poisson distribution
 *
 * @param X a matrix of input values.
 * @param rate_par the rate parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::qpois(X,4);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
qpois(const EigenMat<eT,iTr,iTc>& X, const T1 rate_par)
{
    EIGEN_DIST_FN(qpois_vec,rate_par);
}
#endif
