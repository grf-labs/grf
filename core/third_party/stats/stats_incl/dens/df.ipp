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
 * pdf of the F-distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
df_compute_adj(const T x, const T ratio_abx)
noexcept
{
    return( ratio_abx * (T(1) - x*ratio_abx) );
}

template<typename T>
statslib_constexpr
T
df_log_check(const T x, const T a_par, const T b_par, const T abx, const bool log_form)
noexcept
{
    return( log_form == true ? \
            // if
                dbeta(abx/(T(1)+abx),a_par,b_par,true) \
                + stmath::log(df_compute_adj(x,(a_par/b_par)/(T(1) + abx))) :
            // else
                dbeta(abx/(T(1)+abx),a_par,b_par,false) \
                * df_compute_adj(x,(a_par/b_par)/(T(1) + abx)) );
}

template<typename T>
statslib_constexpr
T
df_limit_vals_dof(const T x, const T df1_par, const T df2_par, const bool log_form)
noexcept
{
    return( // df1 == +Inf and df2 == +Inf
            GCINT::all_posinf(df1_par,df2_par) ? \
                x == T(1) ? \
                    STLIM<T>::infinity() :
                    log_zero_if<T>(log_form) :
            // df1 == +Inf
            GCINT::is_posinf(df1_par) ? \
                // log_form ?
                //     dgamma(1/x,df2_par/2,2/df2_par,true) - 2*stmath::log(x) :
                //     dgamma(1/x,df2_par/2,2/df2_par,false) / (x*x) :
                dinvgamma(x,df2_par/2,df2_par/2,log_form) :
            // df2 == +Inf
                dgamma(x,df1_par/2,2/df1_par,log_form) );
}

template<typename T>
statslib_constexpr
T
df_limit_vals_x(const T df1_par)
noexcept
{
    return( df1_par < T(2) ? \
                STLIM<T>::infinity() :
            //
            df1_par == T(2) ? \
                T(1) :
                T(0) );
}

template<typename T>
statslib_constexpr
T
df_reparam(const T x, const T df1_par, const T df2_par, const bool log_form)
noexcept
{
    return( !f_sanity_check(x,df1_par,df2_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            x < T(0) ? \
                log_zero_if<T>(log_form) :
            x == T(0) ? \
                log_if(df_limit_vals_x(df1_par), log_form) :
            //
            GCINT::any_posinf(df1_par,df2_par) ? \
                df_limit_vals_dof(x,df1_par,df2_par,log_form) :
            //
            df_log_check(x,df1_par/T(2),df2_par/T(2),df1_par*x/df2_par,log_form) );
}

template<typename T1, typename T2, typename T3, typename TC = common_return_t<T1,T2,T3>>
statslib_constexpr
TC
df_type_check(const T1 x, const T2 df1_par, const T3 df2_par, const bool log_form)
noexcept
{
    return df_reparam(static_cast<TC>(x),static_cast<TC>(df1_par),
                      static_cast<TC>(df2_par),log_form);
}

}

/**
 * @brief Density function of the F-distribution
 *
 * @param x a real-valued input.
 * @param df1_par a degrees of freedom parameter, a real-valued input.
 * @param df2_par a degrees of freedom parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return the density function evaluated at \c x.
 * 
 * Example:
 * \code{.cpp} stats::df(1.5,10.0,12.0,false); \endcode
 */

template<typename T1, typename T2, typename T3>
statslib_constexpr
common_return_t<T1,T2,T3>
df(const T1 x, const T2 df1_par, const T3 df2_par, const bool log_form)
noexcept
{
    return internal::df_type_check(x,df1_par,df2_par,log_form);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
void
df_vec(const eT* __stats_pointer_settings__ vals_in, const T1 df1_par, const T2 df2_par, const bool log_form, 
             rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(df,vals_in,vals_out,num_elem,df1_par,df2_par,log_form);
}
#endif

}

/**
 * @brief Density function of the F-distribution
 *
 * @param x a standard vector.
 * @param df1_par a degrees of freedom parameter, a real-valued input.
 * @param df2_par a degrees of freedom parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a vector of density function values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<double> x = {0.3, 0.5, 0.9};
 * stats::df(x,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
std::vector<rT>
df(const std::vector<eT>& x, const T1 df1_par, const T2 df2_par, const bool log_form)
{
    STDVEC_DIST_FN(df_vec,df1_par,df2_par,log_form);
}
#endif

/**
 * @brief Density function of the F-distribution
 *
 * @param X a matrix of input values.
 * @param df1_par a degrees of freedom parameter, a real-valued input.
 * @param df2_par a degrees of freedom parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * arma::mat X = { {0.2,  0.7,  0.1},
 *                 {0.9, -0.3,  1.3} };
 * stats::df(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
ArmaMat<rT>
df(const ArmaMat<eT>& X, const T1 df1_par, const T2 df2_par, const bool log_form)
{
    ARMA_DIST_FN(df_vec,df1_par,df2_par,log_form);
}

template<typename mT, typename tT, typename T1, typename T2>
statslib_inline
mT
df(const ArmaGen<mT,tT>& X, const T1 df1_par, const T2 df2_par, const bool log_form)
{
    return df(X.eval(),df1_par,df2_par,log_form);
}
#endif

/**
 * @brief Density function of the F-distribution
 *
 * @param X a matrix of input values.
 * @param df1_par a degrees of freedom parameter, a real-valued input.
 * @param df2_par a degrees of freedom parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::df(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
df(const BlazeMat<eT,To>& X, const T1 df1_par, const T2 df2_par, const bool log_form)
{
    BLAZE_DIST_FN(df,df1_par,df2_par,log_form);
}
#endif

/**
 * @brief Density function of the F-distribution
 *
 * @param X a matrix of input values.
 * @param df1_par a degrees of freedom parameter, a real-valued input.
 * @param df2_par a degrees of freedom parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::df(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
df(const EigenMat<eT,iTr,iTc>& X, const T1 df1_par, const T2 df2_par, const bool log_form)
{
    EIGEN_DIST_FN(df_vec,df1_par,df2_par,log_form);
}
#endif
