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
 * cdf of the F distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
pf_compute(const T x, const T a_par, const T b_par)
{
    return gcem::incomplete_beta(a_par,b_par, x / (T(1) + x));
}

template<typename T>
statslib_constexpr
T
pf_limit_vals_dof(const T x, const T df1_par, const T df2_par, const bool log_form)
noexcept
{
    return( // df1 == +Inf and df2 == +Inf
            GCINT::all_posinf(df1_par,df2_par) ? \
                x > T(1) ? \
                    log_one_if<T>(log_form) :
                x == T(1) ? \
                    log_if(T(0.5),log_form) :
                    log_zero_if<T>(log_form) :
            // df1 == +Inf
            GCINT::is_posinf(df1_par) ? \
                T(1) - pchisq(df2_par/x,df2_par,log_form) :
            // df2 == +Inf
                pchisq(x*df1_par,df1_par,log_form) );
}

template<typename T>
statslib_constexpr
T
pf_vals_check(const T x, const T df1_par, const T df2_par, const bool log_form)
{
    return( !f_sanity_check(x,df1_par,df2_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            STLIM<T>::epsilon() > x ? \
                log_zero_if<T>(log_form) :
            //
            GCINT::any_posinf(df1_par,df2_par) ? \
                pf_limit_vals_dof(x,df1_par,df2_par,log_form) :
            //
            log_if(pf_compute(df1_par*x/df2_par,df1_par/T(2),df2_par/T(2)), log_form) );
}

template<typename T1, typename T2, typename T3, typename TC = common_return_t<T1,T2,T3>>
statslib_constexpr
TC
pf_type_check(const T1 x, const T2 df1_par, const T3 df2_par, const bool log_form)
noexcept
{
    return pf_vals_check(static_cast<TC>(x),static_cast<TC>(df1_par),
                         static_cast<TC>(df2_par),log_form);
}

}

/**
 * @brief Distribution function of the F-distribution
 *
 * @param x a real-valued input.
 * @param df1_par a degrees of freedom parameter, a real-valued input.
 * @param df2_par a degrees of freedom parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return the cumulative distribution function evaluated at \c x.
 * 
 * Example:
 * \code{.cpp} stats::pf(1.5,10.0,12.0,false); \endcode
 */

template<typename T1, typename T2, typename T3>
statslib_constexpr
common_return_t<T1,T2,T3>
pf(const T1 x, const T2 df1_par, const T3 df2_par, const bool log_form)
noexcept
{
    return internal::pf_type_check(x,df1_par,df2_par,log_form);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
void
pf_vec(const eT* __stats_pointer_settings__ vals_in, const T1 df1_par, const T2 df2_par, const bool log_form, 
             rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(pf,vals_in,vals_out,num_elem,df1_par,df2_par,log_form);
}
#endif

}

/**
 * @brief Distribution function of the Beta distribution
 *
 * @param x a standard vector.
 * @param df1_par a degrees of freedom parameter, a real-valued input.
 * @param df2_par a degrees of freedom parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return a vector of CDF values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<double> x = {0.3, 0.5, 0.9};
 * stats::pf(x,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
std::vector<rT>
pf(const std::vector<eT>& x, const T1 df1_par, const T2 df2_par, const bool log_form)
{
    STDVEC_DIST_FN(pf_vec,df1_par,df2_par,log_form);
}
#endif

/**
 * @brief Distribution function of the Beta distribution
 *
 * @param X a matrix of input values.
 * @param df1_par a degrees of freedom parameter, a real-valued input.
 * @param df2_par a degrees of freedom parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return a matrix of CDF values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * arma::mat X = { {0.2,  0.7,  0.1},
 *                 {0.9, -0.3,  1.3} };
 * stats::pf(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
ArmaMat<rT>
pf(const ArmaMat<eT>& X, const T1 df1_par, const T2 df2_par, const bool log_form)
{
    ARMA_DIST_FN(pf_vec,df1_par,df2_par,log_form);
}

template<typename mT, typename tT, typename T1, typename T2>
statslib_inline
mT
pf(const ArmaGen<mT,tT>& X, const T1 df1_par, const T2 df2_par, const bool log_form)
{
    return pf(X.eval(),df1_par,df2_par,log_form);
}
#endif

/**
 * @brief Distribution function of the Beta distribution
 *
 * @param X a matrix of input values.
 * @param df1_par a degrees of freedom parameter, a real-valued input.
 * @param df2_par a degrees of freedom parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return a matrix of CDF values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::pf(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
pf(const BlazeMat<eT,To>& X, const T1 df1_par, const T2 df2_par, const bool log_form)
{
    BLAZE_DIST_FN(pf_vec,df1_par,df2_par,log_form);
}
#endif

/**
 * @brief Distribution function of the Beta distribution
 *
 * @param X a matrix of input values.
 * @param df1_par a degrees of freedom parameter, a real-valued input.
 * @param df2_par a degrees of freedom parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return a matrix of CDF values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::pf(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
pf(const EigenMat<eT,iTr,iTc>& X, const T1 df1_par, const T2 df2_par, const bool log_form)
{
    EIGEN_DIST_FN(pf_vec,df1_par,df2_par,log_form);
}
#endif
