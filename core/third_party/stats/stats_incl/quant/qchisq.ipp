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
 * quantile function of the chi-squared distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
qchisq_compute(const T p, const T dof_par)
noexcept
{
    return( T(2)*gcem::incomplete_gamma_inv(dof_par/T(2),p) );
}

template<typename T>
statslib_constexpr
T
qchisq_vals_check(const T p, const T dof_par)
noexcept
{
    return( !chisq_sanity_check(dof_par) ? \
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
            qchisq_compute(p,dof_par) );
}

template<typename T1, typename T2, typename TC = common_return_t<T1,T2>>
statslib_constexpr
TC
qchisq_type_check(const T1 p, const T2 dof_par)
noexcept
{
    return qchisq_vals_check(static_cast<TC>(p),static_cast<TC>(dof_par));
}

}

/**
 * @brief Quantile function of the Chi-squared distribution
 *
 * @param p a real-valued input.
 * @param dof_par the degrees of freedom parameter, a real-valued input.
 *
 * @return the quantile function evaluated at \c x.
 * 
 * Example:
 * \code{.cpp} stats::qchisq(0.5,5); \endcode
 */

template<typename T1, typename T2>
statslib_constexpr
common_return_t<T1,T2>
qchisq(const T1 p, const T2 dof_par)
noexcept
{
    return internal::qchisq_type_check(p,dof_par);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename rT>
statslib_inline
void
qchisq_vec(const eT* __stats_pointer_settings__ vals_in, const T1 dof_par, 
                 rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(qchisq,vals_in,vals_out,num_elem,dof_par);
}
#endif

}

/**
 * @brief Quantile function of the Chi-squared distribution
 *
 * @param x a standard vector.
 * @param dof_par the degrees of freedom parameter, a real-valued input.
 *
 * @return a vector of quantile values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<double> x = {0.3, 0.5, 0.8};
 * stats::qchisq(x,4);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
std::vector<rT>
qchisq(const std::vector<eT>& x, const T1 dof_par)
{
    STDVEC_DIST_FN(qchisq_vec,dof_par);
}
#endif

/**
 * @brief Quantile function of the Chi-squared distribution
 *
 * @param X a matrix of input values.
 * @param dof_par the degrees of freedom parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * arma::mat X = { {0.2, 0.7, 0.9},
 *                 {0.1, 0.8, 0.3} };
 * stats::qchisq(X,4);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
ArmaMat<rT>
qchisq(const ArmaMat<eT>& X, const T1 dof_par)
{
    ARMA_DIST_FN(qchisq_vec,dof_par);
}

template<typename mT, typename tT, typename T1>
statslib_inline
mT
qchisq(const ArmaGen<mT,tT>& X, const T1 dof_par)
{
    return qchisq(X.eval(),dof_par);
}
#endif

/**
 * @brief Quantile function of the Chi-squared distribution
 *
 * @param X a matrix of input values.
 * @param dof_par the degrees of freedom parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::qchisq(X,4);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
qchisq(const BlazeMat<eT,To>& X, const T1 dof_par)
{
    BLAZE_DIST_FN(qchisq_vec,dof_par);
}
#endif

/**
 * @brief Quantile function of the Chi-squared distribution
 *
 * @param X a matrix of input values.
 * @param dof_par the degrees of freedom parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::qchisq(X,4);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
qchisq(const EigenMat<eT,iTr,iTc>& X, const T1 dof_par)
{
    EIGEN_DIST_FN(qchisq_vec,dof_par);
}
#endif
