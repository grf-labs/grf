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
 * quantile function of the F-distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
qf_compute_adj(const T I_inv_val, const T ab_ratio)
noexcept
{
    return( I_inv_val / (ab_ratio*(T(1) - I_inv_val)) );
}

template<typename T>
statslib_constexpr
T
qf_compute(const T p, const T a_par, const T b_par)
noexcept
{
    return qf_compute_adj(gcem::incomplete_beta_inv(a_par,b_par,p),a_par/b_par);
}

template<typename T>
statslib_constexpr
T
qf_limit_vals_dof(const T p, const T df1_par, const T df2_par)
noexcept
{
    return( // df1 == +Inf and df2 == +Inf
            GCINT::all_posinf(df1_par,df2_par) ? \
                T(1) :
            // df1 == +Inf
            GCINT::is_posinf(df1_par) ? \
                df2_par / qchisq(T(1)-p,df2_par) :
            // df2 == +Inf
                qchisq(p,df1_par) / df1_par );
}

template<typename T>
statslib_constexpr
T
qf_vals_check(const T p, const T df1_par, const T df2_par)
noexcept
{
    return( !f_sanity_check(df1_par,df2_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            !prob_val_check(p) ? \
                STLIM<T>::quiet_NaN() :
            //
            p == T(0) ? \
                T(0) :
            p == T(1) ? \
                STLIM<T>::infinity() :
            // 0 < p < 1
            GCINT::any_posinf(df1_par,df2_par) ? \
                qf_limit_vals_dof(p,df1_par,df2_par) :
            //
            qf_compute(p,df1_par/T(2),df2_par/T(2)) );
}

template<typename T1, typename T2, typename T3, typename TC = common_return_t<T1,T2,T3>>
statslib_constexpr
TC
qf_type_check(const T1 p, const T2 df1_par, const T3 df2_par)
noexcept
{
    return qf_vals_check(static_cast<TC>(p),static_cast<TC>(df1_par),static_cast<TC>(df2_par));
}

}

/**
 * @brief Quantile function of the F-distribution
 *
 * @param p a real-valued input.
 * @param df1_par a degrees of freedom parameter, a real-valued input.
 * @param df2_par a degrees of freedom parameter, a real-valued input.
 *
 * @return the quantile function evaluated at \c p.
 * 
 * Example:
 * \code{.cpp} stats::qf(0.5,10.0,12.0); \endcode
 */

template<typename T1, typename T2, typename T3>
statslib_constexpr
common_return_t<T1,T2,T3>
qf(const T1 p, const T2 df1_par, const T3 df2_par)
noexcept
{
    return internal::qf_type_check(p,df1_par,df2_par);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
void
qf_vec(const eT* __stats_pointer_settings__ vals_in, const T1 df1_par, const T2 df2_par, 
             rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(qf,vals_in,vals_out,num_elem,df1_par,df2_par);
}
#endif

}

/**
 * @brief Quantile function of the F-distribution
 *
 * @param x a standard vector.
 * @param df1_par a degrees of freedom parameter, a real-valued input.
 * @param df2_par a degrees of freedom parameter, a real-valued input.
 *
 * @return a vector of quantile values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<double> x = {0.3, 0.5, 0.9};
 * stats::qf(x,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
std::vector<rT>
qf(const std::vector<eT>& x, const T1 df1_par, const T2 df2_par)
{
    STDVEC_DIST_FN(qf_vec,df1_par,df2_par);
}
#endif

/**
 * @brief Quantile function of the F-distribution
 *
 * @param X a matrix of input values.
 * @param df1_par a degrees of freedom parameter, a real-valued input.
 * @param df2_par a degrees of freedom parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * arma::mat X = { {0.2,  0.7,  0.1},
 *                 {0.9,  0.3,  0.87} };
 * stats::qf(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
ArmaMat<rT>
qf(const ArmaMat<eT>& X, const T1 df1_par, const T2 df2_par)
{
    ARMA_DIST_FN(qf_vec,df1_par,df2_par);
}

template<typename mT, typename tT, typename T1, typename T2>
statslib_inline
mT
qf(const ArmaGen<mT,tT>& X, const T1 df1_par, const T2 df2_par)
{
    return qf(X.eval(),df1_par,df2_par);
}
#endif

/**
 * @brief Quantile function of the F-distribution
 *
 * @param X a matrix of input values.
 * @param df1_par a degrees of freedom parameter, a real-valued input.
 * @param df2_par a degrees of freedom parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::qf(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
qf(const BlazeMat<eT,To>& X, const T1 df1_par, const T2 df2_par)
{
    BLAZE_DIST_FN(qf_vec,df1_par,df2_par);
}
#endif

/**
 * @brief Quantile function of the F-distribution
 *
 * @param X a matrix of input values.
 * @param df1_par a degrees of freedom parameter, a real-valued input.
 * @param df2_par a degrees of freedom parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::qf(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
qf(const EigenMat<eT,iTr,iTc>& X, const T1 df1_par, const T2 df2_par)
{
    EIGEN_DIST_FN(qf_vec,df1_par,df2_par);
}
#endif
