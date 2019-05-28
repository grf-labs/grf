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
 * pdf of the Beta distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
dbeta_log_compute(const T x, const T a_par, const T b_par)
noexcept
{
    return( - (stmath::lgamma(a_par) + stmath::lgamma(b_par) - stmath::lgamma(a_par+b_par)) \
            + (a_par - T(1))*stmath::log(x) + (b_par - T(1))*stmath::log(T(1) - x) );
}

template<typename T>
statslib_constexpr
T
dbeta_limit_vals(const T x, const T a_par, const T b_par)
noexcept
{
    return( a_par == T(0) && b_par == T(0) ? \
                x == T(0) || x == T(1) ? \
                    STLIM<T>::infinity() :
                    T(0) : 
            //
            a_par == T(0) || (gcem::internal::is_posinf(b_par) && !gcem::internal::is_posinf(a_par)) ? \
                x == T(0) ? \
                    STLIM<T>::infinity() :
                    T(0) :
            //
            b_par == T(0) || (gcem::internal::is_posinf(a_par) && !gcem::internal::is_posinf(b_par)) ? \
                x == T(1) ? \
                    STLIM<T>::infinity() :
                    T(0) :
            //
            x == T(0.5) ? \
                STLIM<T>::infinity() : 
                T(0) );
}

template<typename T>
statslib_constexpr
T
dbeta_vals_check(const T x, const T a_par, const T b_par, const bool log_form)
noexcept
{
    return( !beta_sanity_check(x,a_par,b_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            x < T(0) || x > T(1) ? \
                log_zero_if<T>(log_form) :
            //
            (a_par == T(0) || b_par == T(0) || gcem::internal::any_posinf(a_par,b_par)) ? \
                log_if(dbeta_limit_vals(x,a_par,b_par),log_form) :
            //
            x == T(0) ? \
                a_par < T(1) ? \
                    STLIM<T>::infinity() :
                a_par > T(1) ? \
                    log_zero_if<T>(log_form) :
                    log_if(b_par,log_form) :
            x == T(1) ? \
                b_par < T(1) ? \
                    STLIM<T>::infinity() :
                b_par > T(1) ? \
                    log_zero_if<T>(log_form) :
                    log_if(a_par,log_form) :
            //
            exp_if(dbeta_log_compute(x,a_par,b_par), !log_form) );
}

template<typename T1, typename T2, typename T3, typename TC = common_return_t<T1,T2,T3>>
statslib_constexpr
TC
dbeta_type_check(const T1 x, const T2 a_par, const T3 b_par, const bool log_form)
noexcept
{
    return dbeta_vals_check(static_cast<TC>(x),static_cast<TC>(a_par),
                            static_cast<TC>(b_par),log_form);
}

}

/**
 * @brief Density function of the Beta distribution
 *
 * @param x a real-valued input.
 * @param a_par a real-valued shape parameter.
 * @param b_par a real-valued shape parameter.
 * @param log_form return the log-density or the true form.
 *
 * @return the density function evaluated at \c x.
 * 
 * Example:
 * \code{.cpp} stats::dbeta(0.5,3.0,2.0,false); \endcode
 */

template<typename T1, typename T2, typename T3>
statslib_constexpr
common_return_t<T1,T2,T3>
dbeta(const T1 x, const T2 a_par, const T3 b_par, const bool log_form)
noexcept
{
    return internal::dbeta_type_check(x,a_par,b_par,log_form);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
void
dbeta_vec(const eT* __stats_pointer_settings__ vals_in, const T1 a_par, const T2 b_par, const bool log_form, 
                rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(dbeta,vals_in,vals_out,num_elem,a_par,b_par,log_form);
}
#endif

}

/**
 * @brief Density function of the Beta distribution
 *
 * @param x a standard vector.
 * @param a_par a real-valued shape parameter.
 * @param b_par a real-valued shape parameter.
 * @param log_form return the log-density or the true form.
 *
 * @return a vector of density function values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<double> x = {0.3, 0.5, 0.9};
 * stats::dbeta(x,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
std::vector<rT>
dbeta(const std::vector<eT>& x, const T1 a_par, const T2 b_par, const bool log_form)
{
    STDVEC_DIST_FN(dbeta_vec,a_par,b_par,log_form);
}
#endif

/**
 * @brief Density function of the Beta distribution
 *
 * @param X a matrix of input values.
 * @param a_par a real-valued shape parameter.
 * @param b_par a real-valued shape parameter.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * arma::mat X = { {0.2,  0.7,  0.1},
 *                 {0.9,  0.3,  0.87} };
 * stats::dbeta(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
ArmaMat<rT>
dbeta(const ArmaMat<eT>& X, const T1 a_par, const T2 b_par, const bool log_form)
{
    ARMA_DIST_FN(dbeta_vec,a_par,b_par,log_form);
}

template<typename mT, typename tT, typename T1, typename T2>
statslib_inline
mT
dbeta(const ArmaGen<mT,tT>& X, const T1 a_par, const T2 b_par, const bool log_form)
{
    return dbeta(X.eval(),a_par,b_par,log_form);
}
#endif

/**
 * @brief Density function of the Beta distribution
 *
 * @param X a matrix of input values.
 * @param a_par a real-valued shape parameter.
 * @param b_par a real-valued shape parameter.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::dbeta(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
dbeta(const BlazeMat<eT,To>& X, const T1 a_par, const T2 b_par, const bool log_form)
{
    BLAZE_DIST_FN(dbeta,a_par,b_par,log_form);
}
#endif

/**
 * @brief Density function of the Beta distribution
 *
 * @param X a matrix of input values.
 * @param a_par a real-valued shape parameter.
 * @param b_par a real-valued shape parameter.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::dbeta(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
dbeta(const EigenMat<eT,iTr,iTc>& X, const T1 a_par, const T2 b_par, const bool log_form)
{
    EIGEN_DIST_FN(dbeta_vec,a_par,b_par,log_form);
}
#endif
