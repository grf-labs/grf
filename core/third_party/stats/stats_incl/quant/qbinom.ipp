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
 * quantile function of the Binomial distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
qbinom_recur(const T p, const llint_t n_trials_par, const T prob_par, const T value, const llint_t count)
noexcept
{
    return( value < p ? \
                qbinom_recur(p,n_trials_par,prob_par, pbinom(count,n_trials_par,prob_par), count + 1) :
                count - llint_t(1) );
}

template<typename T>
statslib_constexpr
T
qbinom_vals_check(const T p, const llint_t n_trials_par, const T prob_par)
noexcept
{
    return( !binom_sanity_check(n_trials_par,prob_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            !prob_val_check(p) ? \
                STLIM<T>::quiet_NaN() :
            //
            p == T(0) ? \
                T(0) :
            p == T(1) ? \
                static_cast<T>(n_trials_par) :
            //
            qbinom_recur(p,n_trials_par,prob_par,T(0),llint_t(0)) );
}

template<typename T1, typename T2, typename TC = common_return_t<T1,T2>>
statslib_constexpr
TC
qbinom_type_check(const T1 p, const llint_t n_trials_par, const T2 prob_par)
noexcept
{
    return qbinom_vals_check(static_cast<TC>(p),n_trials_par,static_cast<TC>(prob_par));
}

}

/**
 * @brief Quantile function of the Binomial distribution
 *
 * @param p a real-valued input.
 * @param n_trials_par the number of trials, a non-negative integral-valued input.
 * @param prob_par the probability parameter, a real-valued input.
 *
 * @return the quantile function evaluated at \c p.
 * 
 * Example:
 * \code{.cpp} stats::qbinom(0.4,4,0.4); \endcode
 */

template<typename T1, typename T2>
statslib_constexpr
common_return_t<T1,T2>
qbinom(const T1 p, const llint_t n_trials_par, const T2 prob_par)
noexcept
{
    return internal::qbinom_type_check(p,n_trials_par,prob_par);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename rT>
statslib_inline
void
qbinom_vec(const eT* __stats_pointer_settings__ vals_in, const llint_t n_trials_par, const T1 prob_par, 
                 rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(qbinom,vals_in,vals_out,num_elem,n_trials_par,prob_par);
}
#endif

}

/**
 * @brief Quantile function of the Binomial distribution
 *
 * @param x a standard vector.
 * @param n_trials_par the number of trials, a non-negative integral-valued input.
 * @param prob_par the probability parameter, a real-valued input.
 *
 * @return a vector of quantile values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<int> x = {2, 3, 4};
 * stats::qbinom(x,5,0.5);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
std::vector<rT>
qbinom(const std::vector<eT>& x, const llint_t n_trials_par, const T1 prob_par)
{
    STDVEC_DIST_FN(qbinom_vec,n_trials_par,prob_par);
}
#endif

/**
 * @brief Quantile function of the Binomial distribution
 *
 * @param X a matrix of input values.
 * @param n_trials_par the number of trials, a non-negative integral-valued input.
 * @param prob_par the probability parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::qbinom(X,5,0.5);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
ArmaMat<rT>
qbinom(const ArmaMat<eT>& X, const llint_t n_trials_par, const T1 prob_par)
{
    ARMA_DIST_FN(qbinom_vec,n_trials_par,prob_par);
}

template<typename mT, typename tT, typename T1>
statslib_inline
mT
qbinom(const ArmaGen<mT,tT>& X, const llint_t n_trials_par, const T1 prob_par)
{
    return qbinom(X.eval(),n_trials_par,prob_par);
}
#endif

/**
 * @brief Quantile function of the Binomial distribution
 *
 * @param X a matrix of input values.
 * @param n_trials_par the number of trials, a non-negative integral-valued input.
 * @param prob_par the probability parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::qbinom(X,5,0.5);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
qbinom(const BlazeMat<eT,To>& X, const llint_t n_trials_par, const T1 prob_par)
{
    BLAZE_DIST_FN(qbinom_vec,n_trials_par,prob_par);
}
#endif

/**
 * @brief Quantile function of the Binomial distribution
 *
 * @param X a matrix of input values.
 * @param n_trials_par the number of trials, a non-negative integral-valued input.
 * @param prob_par the probability parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::qbinom(X,5,0.5);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
qbinom(const EigenMat<eT,iTr,iTc>& X, const llint_t n_trials_par, const T1 prob_par)
{
    EIGEN_DIST_FN(qbinom_vec,n_trials_par,prob_par);
}
#endif
