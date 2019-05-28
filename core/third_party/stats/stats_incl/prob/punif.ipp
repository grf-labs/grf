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
 * cdf of the uniform distribution
 */

//
// single input

namespace internal
{

template<typename T>
statslib_constexpr
T
punif_compute(const T x, const T a_par, const T b_par)
{
    return( x <= a_par ? \
                T(0) :
            x >= b_par ? \
                T(1) :
            //
                (x-a_par) / (b_par-a_par) );
}

template<typename T>
statslib_constexpr
T
punif_vals_check(const T x, const T a_par, const T b_par, const bool log_form)
{
    return( !unif_sanity_check(x,a_par,b_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            GCINT::all_inf(a_par,b_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            log_if(punif_compute(x,a_par,b_par), log_form) );
}

template<typename T1, typename T2, typename T3, typename TC = common_return_t<T1,T2,T3>>
statslib_constexpr
TC
punif_type_check(const T1 x, const T2 a_par, const T3 b_par, const bool log_form)
noexcept
{
    return punif_vals_check(static_cast<TC>(x),static_cast<TC>(a_par),
                           static_cast<TC>(b_par),log_form);
}

}

/**
 * @brief Distribution function of the Uniform distribution
 *
 * @param x a real-valued input.
 * @param a_par the lower bound parameter, a real-valued input.
 * @param b_par the upper bound parameter, a real-valued input.
 * @param log_form return the log-probability or the true form.
 *
 * @return the cumulative distribution function evaluated at \c x.
 * 
 * Example:
 * \code{.cpp} stats::punif(0.5,-1.0,2.0,false); \endcode
 */

template<typename T1, typename T2, typename T3>
statslib_constexpr
common_return_t<T1,T2,T3>
punif(const T1 x, const T2 a_par, const T3 b_par, const bool log_form)
noexcept
{
    return internal::punif_type_check(x,a_par,b_par,log_form);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
void
punif_vec(const eT* __stats_pointer_settings__ vals_in, const T1 a_par, const T2 b_par, const bool log_form, 
                rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(punif,vals_in,vals_out,num_elem,a_par,b_par,log_form);
}
#endif

}

/**
 * @brief Distribution function of the Uniform distribution
 *
 * @param x a standard vector.
 * @param a_par a real-valued shape parameter.
 * @param b_par a real-valued shape parameter.
 * @param log_form return the log-probability or the true form.
 *
 * @return a vector of CDF values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<double> x = {0.3, 0.5, 0.9};
 * stats::punif(x,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
std::vector<rT>
punif(const std::vector<eT>& x, const T1 a_par, const T2 b_par, const bool log_form)
{
    STDVEC_DIST_FN(punif_vec,a_par,b_par,log_form);
}
#endif

/**
 * @brief Distribution function of the Uniform distribution
 *
 * @param X a matrix of input values.
 * @param a_par a real-valued shape parameter.
 * @param b_par a real-valued shape parameter.
 * @param log_form return the log-probability or the true form.
 *
 * @return a matrix of CDF values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * arma::mat X = { {0.2,  0.7,  0.1},
 *                 {0.9, -0.3,  1.3} };
 * stats::punif(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
ArmaMat<rT>
punif(const ArmaMat<eT>& X, const T1 a_par, const T2 b_par, const bool log_form)
{
    ARMA_DIST_FN(punif_vec,a_par,b_par,log_form);
}

template<typename mT, typename tT, typename T1, typename T2>
statslib_inline
mT
punif(const ArmaGen<mT,tT>& X, const T1 a_par, const T2 b_par, const bool log_form)
{
    return punif(X.eval(),a_par,b_par,log_form);
}
#endif

/**
 * @brief Distribution function of the Uniform distribution
 *
 * @param X a matrix of input values.
 * @param a_par a real-valued shape parameter.
 * @param b_par a real-valued shape parameter.
 * @param log_form return the log-probability or the true form.
 *
 * @return a matrix of CDF values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::punif(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
punif(const BlazeMat<eT,To>& X, const T1 a_par, const T2 b_par, const bool log_form)
{
    BLAZE_DIST_FN(punif_vec,a_par,b_par,log_form);
}
#endif

/**
 * @brief Distribution function of the Uniform distribution
 *
 * @param X a matrix of input values.
 * @param a_par a real-valued shape parameter.
 * @param b_par a real-valued shape parameter.
 * @param log_form return the log-probability or the true form.
 *
 * @return a matrix of CDF values corresponding to the elements of \c X.
 *
 * Example:
 * \code{.cpp}
 * stats::punif(X,3.0,2.0,false);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
punif(const EigenMat<eT,iTr,iTc>& X, const T1 a_par, const T2 b_par, const bool log_form)
{
    EIGEN_DIST_FN(punif_vec,a_par,b_par,log_form);
}
#endif
