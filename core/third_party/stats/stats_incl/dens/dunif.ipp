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
 * pdf of the uniform distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
dunif_log_check(const T a_par, const T b_par, const bool log_form)
noexcept
{
    return( log_form == true ? \
                - stmath::log(b_par - a_par) : 
                T(1) / (b_par - a_par) );
}

template<typename T>
statslib_constexpr
T
dunif_vals_check(const T x, const T a_par, const T b_par, const bool log_form)
noexcept
{
    return( !unif_sanity_check(x,a_par,b_par) ? \
                STLIM<T>::quiet_NaN() :
            // cases: a_par == -Inf, b_par finite; a_par finite, b_par == +Inf; 
            // and a_par == -Inf, b_par == +Inf
            GCINT::any_inf(a_par,b_par) ? \
                T(0) :
            //
            x < a_par || x > b_par ? \
                log_zero_if<T>(log_form) :
            //
            dunif_log_check(a_par,b_par,log_form) );
}

template<typename T1, typename T2, typename T3, typename TC = common_return_t<T1,T2,T3>>
statslib_constexpr
TC
dunif_type_check(const T1 x, const T2 a_par, const T3 b_par, const bool log_form)
noexcept
{
    return dunif_vals_check(static_cast<TC>(x),static_cast<TC>(a_par),
                            static_cast<TC>(b_par),log_form);
}

}

/**
 * @brief Density function of the Uniform distribution
 *
 * @param x a real-valued input.
 * @param a_par the lower bound parameter, a real-valued input.
 * @param b_par the upper bound parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return the density function evaluated at \c x.
 * 
 * Example:
 * \code{.cpp} stats::dunif(0.5,-1.0,2.0,false); \endcode
 */

template<typename T1, typename T2, typename T3>
statslib_constexpr
common_return_t<T1,T2,T3>
dunif(const T1 x, const T2 a_par, const T3 b_par, const bool log_form)
noexcept
{
    return internal::dunif_type_check(x,a_par,b_par,log_form);
}

template<typename T>
statslib_constexpr
return_t<T>
dunif(const T x)
noexcept
{
    return dunif(x,T(0),T(1));
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
void
dunif_vec(const eT* __stats_pointer_settings__ vals_in, const T1 a_par, const T2 b_par, const bool log_form, 
                rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(dunif,vals_in,vals_out,num_elem,a_par,b_par,log_form);
}
#endif

}

/**
 * @brief Density function of the Uniform distribution
 *
 * @param x a standard vector.
 * @param a_par the lower bound parameter, a real-valued input.
 * @param b_par the upper bound parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a vector of density function values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<double> x = {-2.0, 0.0, 2.0};
 * stats::dunif(x,-1.0,3.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
std::vector<rT>
dunif(const std::vector<eT>& x, const T1 a_par, const T2 b_par, const bool log_form)
{
    STDVEC_DIST_FN(dunif_vec,a_par,b_par,log_form);
}
#endif

/**
 * @brief Density function of the Uniform distribution
 *
 * @param X a matrix of input values.
 * @param a_par the lower bound parameter, a real-valued input.
 * @param b_par the upper bound parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * arma::mat X = { {0.2,  0.7,  0.1},
 *                 {0.9, -0.3,  1.3} };
 * stats::dunif(X,-1.0,3.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
ArmaMat<rT>
dunif(const ArmaMat<eT>& X, const T1 a_par, const T2 b_par, const bool log_form)
{
    ARMA_DIST_FN(dunif_vec,a_par,b_par,log_form);
}

template<typename mT, typename tT, typename T1, typename T2>
statslib_inline
mT
dunif(const ArmaGen<mT,tT>& X, const T1 a_par, const T2 b_par, const bool log_form)
{
    return dunif(X.eval(),a_par,b_par,log_form);
}
#endif

/**
 * @brief Density function of the Uniform distribution
 *
 * @param X a matrix of input values.
 * @param a_par the lower bound parameter, a real-valued input.
 * @param b_par the upper bound parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::dunif(X,-1.0,3.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
dunif(const BlazeMat<eT,To>& X, const T1 a_par, const T2 b_par, const bool log_form)
{
    BLAZE_DIST_FN(dunif,a_par,b_par,log_form);
}
#endif

/**
 * @brief Density function of the Uniform distribution
 *
 * @param X a matrix of input values.
 * @param a_par the lower bound parameter, a real-valued input.
 * @param b_par the upper bound parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::dunif(X,-1.0,3.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
dunif(const EigenMat<eT,iTr,iTc>& X, const T1 a_par, const T2 b_par, const bool log_form)
{
    EIGEN_DIST_FN(dunif_vec,a_par,b_par,log_form);
}
#endif
