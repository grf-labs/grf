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
 * Sample from an exponential distribution
 */

//
// scalar output

namespace internal
{

template<typename T>
statslib_inline
T
rexp_compute(const T rate_par, rand_engine_t& engine)
{
    return( !exp_sanity_check(rate_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            qexp(runif(T(0),T(1),engine),rate_par) );
}

}

/**
 * @brief Random sampling function for the Exponential distribution
 *
 * @param rate_par the rate parameter, a real-valued input.
 * @param engine a random engine, passed by reference.
 *
 * @return a pseudo-random draw from the Exponential distribution.
 *
 * Example:
 * \code{.cpp}
 * stats::rand_engine_t engine(1776);
 * stats::rexp(4,engine);
 * \endcode
 */

template<typename T>
statslib_inline
return_t<T>
rexp(const T rate_par, rand_engine_t& engine)
{
    return internal::rexp_compute(static_cast<return_t<T>>(rate_par),engine);
}

/**
 * @brief Random sampling function for the Exponential distribution
 *
 * @param rate_par the rate parameter, a real-valued input.
 * @param seed_val initialize the random engine with a non-negative integral-valued seed.
 *
 * @return a pseudo-random draw from the Exponential distribution.
 *
 * Example:
 * \code{.cpp}
 * stats::rexp(4,1776);
 * \endcode
 */

template<typename T>
statslib_inline
return_t<T>
rexp(const T rate_par, const ullint_t seed_val)
{
    rand_engine_t engine(seed_val);
    return rexp(rate_par,engine);
}

//
// vector/matrix output

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename T1, typename rT>
statslib_inline
void
rexp_vec(const T1 rate_par, rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    RAND_DIST_FN_VEC(rexp,vals_out,num_elem,rate_par);
}
#endif

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1>
statslib_inline
void
rexp_mat_check(std::vector<eT>& X, const T1 dof_par)
{
    STDVEC_RAND_DIST_FN(rexp,dof_par);
}
#endif

#ifdef STATS_ENABLE_MATRIX_FEATURES
template<typename mT, typename T1>
statslib_inline
void
rexp_mat_check(mT& X, const T1 dof_par)
{
    MAIN_MAT_RAND_DIST_FN(rexp,dof_par);
}
#endif

}

/**
 * @brief Random matrix sampling function for the Exponential distribution
 *
 * @param n the number of output rows
 * @param k the number of output columns
 * @param rate_par the rate parameter, a real-valued input.
 *
 * @return a matrix of pseudo-random draws from the Exponential distribution.
 *
 * Example:
 * \code{.cpp}
 * // std::vector
 * stats::rexp<std::vector<double>>(5,4,4);
 * // Armadillo matrix
 * stats::rexp<arma::mat>(5,4,4);
 * // Blaze dynamic matrix
 * stats::rexp<blaze::DynamicMatrix<double,blaze::columnMajor>>(5,4,4);
 * // Eigen dynamic matrix
 * stats::rexp<Eigen::MatrixXd>(5,4,4);
 * \endcode
 *
 * @note This function requires template instantiation; acceptable output types include: <tt>std::vector</tt> with primitive types (e.g., \c float, \c double, etc.), as well as Armadillo, Blaze, and Eigen dense matrices.
 */

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename mT, typename T1>
statslib_inline
mT
rexp(const ullint_t n, const ullint_t k, const T1 rate_par)
{
    GEN_MAT_RAND_FN(rexp_mat_check,rate_par);
}
#endif
