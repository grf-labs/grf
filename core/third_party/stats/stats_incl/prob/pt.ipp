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
 * cdf of the univariate t distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
pt_compute_main_2(const T z, const T r_par)
{
    return( pbeta(T(1)/z, r_par/T(2), T(0.5)) / T(2) );
}

template<typename T>
statslib_constexpr
T
pt_compute_main_1(const T z, const T r_par)
{
    return( T(0.5) - pbeta(z/(r_par+z), T(0.5), r_par/T(2)) / T(2) );
}

template<typename T>
statslib_constexpr
T
pt_compute_main(const T x, const T r_par)
{
    return( r_par > x*x ? \
                (x > T(0) ? T(1) - pt_compute_main_1(x*x,r_par) : pt_compute_main_1(x*x,r_par)) : 
                (x > T(0) ? T(1) - pt_compute_main_2(T(1) + (x/r_par)*x,r_par) : 
                            pt_compute_main_2(T(1) + (x/r_par)*x,r_par)) );
}

template<typename T>
statslib_constexpr
T
pt_compute(const T x, const T r_par)
{
    return( r_par == T(1) ? \
                pcauchy_compute(x) :
            r_par == T(2) ? \
                T(0.5) + x / (T(2) * stmath::sqrt(x*x + T(2)) ) :
            //
            pt_compute_main(x,T(r_par)) );
}

template<typename T>
statslib_constexpr
T
pt_vals_check(const T x, const T dof_par, const bool log_form)
{
    return( !t_sanity_check(x,dof_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            GCINT::is_posinf(x) ? \
                log_one_if<T>(log_form) :
            GCINT::is_neginf(x) ? \
                log_zero_if<T>(log_form) :
            //
            GCINT::is_posinf(dof_par) ? \
                pnorm(x,T(0),T(1),log_form) :
            //
            log_if(pt_compute(x,dof_par), log_form) );
}

template<typename T1, typename T2, typename TC = common_return_t<T1,T2>>
statslib_constexpr
TC
pt_type_check(const T1 x, const T2 dof_par, const bool log_form)
noexcept
{
    return pt_vals_check(static_cast<TC>(x),static_cast<TC>(dof_par),log_form);
}

}

/**
 * @brief Distribution function of the t-distribution
 *
 * @param x a real-valued input.
 * @param dof_par the degrees of freedom parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return the cumulative distribution function evaluated at \c x.
 * 
 * Example:
 * \code{.cpp} stats::pt(0.37,11,false); \endcode
 */

template<typename T1, typename T2>
statslib_constexpr
common_return_t<T1,T2>
pt(const T1 x, const T2 dof_par, const bool log_form)
noexcept
{
    return internal::pt_type_check(x,dof_par,log_form);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename rT>
statslib_inline
void
pt_vec(const eT* __stats_pointer_settings__ vals_in, const T1 dof_par, const bool log_form, 
             rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(pt,vals_in,vals_out,num_elem,dof_par,log_form);
}
#endif

}

/**
 * @brief Distribution function of the t-distribution
 *
 * @param x a standard vector.
 * @param dof_par the degrees of freedom parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return a vector of CDF values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<double> x = {0.0, 1.0, 2.0};
 * stats::pt(x,4,false);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
std::vector<rT>
pt(const std::vector<eT>& x, const T1 dof_par, const bool log_form)
{
    STDVEC_DIST_FN(pt_vec,dof_par,log_form);
}
#endif

/**
 * @brief Distribution function of the t-distribution
 *
 * @param X a matrix of input values.
 * @param dof_par the degrees of freedom parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return a matrix of CDF values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * arma::mat X = { {0.2, -1.7,  0.1},
 *                 {0.9,  4.0, -0.3} };
 * stats::pt(X,4,false);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
ArmaMat<rT>
pt(const ArmaMat<eT>& X, const T1 dof_par, const bool log_form)
{
    ARMA_DIST_FN(pt_vec,dof_par,log_form);
}

template<typename mT, typename tT, typename T1>
statslib_inline
mT
pt(const ArmaGen<mT,tT>& X, const T1 dof_par, const bool log_form)
{
    return pt(X.eval(),dof_par,log_form);
}
#endif

/**
 * @brief Distribution function of the t-distribution
 *
 * @param X a matrix of input values.
 * @param dof_par the degrees of freedom parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return a matrix of CDF values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::pt(X,4,false);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
pt(const BlazeMat<eT,To>& X, const T1 dof_par, const bool log_form)
{
    BLAZE_DIST_FN(pt_vec,dof_par,log_form);
}
#endif

/**
 * @brief Distribution function of the t-distribution
 *
 * @param X a matrix of input values.
 * @param dof_par the degrees of freedom parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return a matrix of CDF values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::pt(X,4,false);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
pt(const EigenMat<eT,iTr,iTc>& X, const T1 dof_par, const bool log_form)
{
    EIGEN_DIST_FN(pt_vec,dof_par,log_form);
}
#endif
