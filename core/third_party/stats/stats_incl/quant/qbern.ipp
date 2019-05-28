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
 * quantile function of the univariate Bernoulli distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
qbern_compute(const T p, const T prob_par)
noexcept
{
    return( !bern_sanity_check(prob_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            !prob_val_check(p) ? \
                STLIM<T>::quiet_NaN() :
            //
            p > (T(1) - prob_par) ? \
                llint_t(1) : 
                llint_t(0) );
}

template<typename T1, typename T2, typename TC = common_return_t<T1,T2>>
statslib_constexpr
TC
qbern_type_check(const T1 p, const T2 prob_par)
noexcept
{
    return qbern_compute(static_cast<TC>(p),static_cast<TC>(prob_par));
}

}

/**
 * @brief Quantile function of the Bernoulli distribution
 *
 * @param p a real-valued input.
 * @param prob_par the probability parameter, a real-valued input.
 *
 * @return the quantile function evaluated at \c p.
 * 
 * Example:
 * \code{.cpp} stats::qbern(0.5,0.4); \endcode
 */

template<typename T1, typename T2>
statslib_constexpr
common_return_t<T1,T2> // not llint_t so we can return NaN
qbern(const T1 p, const T2 prob_par)
noexcept
{
    return internal::qbern_type_check(p,prob_par);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename rT>
statslib_inline
void
qbern_vec(const eT* __stats_pointer_settings__ vals_in, const T1 prob_par,
                rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(qbern,vals_in,vals_out,num_elem,prob_par);
}
#endif

}

/**
 * @brief Quantile function of the Bernoulli distribution
 *
 * @param x a standard vector.
 * @param prob_par the probability parameter, a real-valued input.
 *
 * @return a vector of quantile values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<int> x = {0.4, 0.5, 0.9};
 * stats::qbern(x,0.5);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
std::vector<rT>
qbern(const std::vector<eT>& x, const T1 prob_par)
{
    STDVEC_DIST_FN(qbern_vec,prob_par);
}
#endif

/**
 * @brief Quantile function of the Bernoulli distribution
 *
 * @param X a matrix of input values.
 * @param prob_par the probability parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * arma::mat X = { {0.4, 0.5, 0.9},
 *                 {0.3, 0.6, 0.7} };
 * stats::qbern(X,0.5);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
ArmaMat<rT>
qbern(const ArmaMat<eT>& X, const T1 prob_par)
{
    ARMA_DIST_FN(qbern_vec,prob_par);
}

template<typename mT, typename tT, typename T1>
statslib_inline
mT
qbern(const ArmaGen<mT,tT>& X, const T1 prob_par)
{
    return qbern(X.eval(),prob_par);
}
#endif

/**
 * @brief Quantile function of the Bernoulli distribution
 *
 * @param X a matrix of input values.
 * @param prob_par the probability parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::qbern(X,0.5);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
qbern(const BlazeMat<eT,To>& X, const T1 prob_par)
{
    BLAZE_DIST_FN(qbern_vec,prob_par);
}
#endif

/**
 * @brief Quantile function of the Bernoulli distribution
 *
 * @param X a matrix of input values.
 * @param prob_par the probability parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::qbern(X,0.5);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
qbern(const EigenMat<eT,iTr,iTc>& X, const T1 prob_par)
{
    EIGEN_DIST_FN(qbern_vec,prob_par);
}
#endif
