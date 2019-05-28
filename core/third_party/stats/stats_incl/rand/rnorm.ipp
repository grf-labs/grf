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
 * Sample from a normal distribution
 */

//
// single draw

namespace internal
{

template<typename T>
statslib_inline
T
rnorm_compute(const T mu_par, const T sigma_par, rand_engine_t& engine)
{
    if (!norm_sanity_check(mu_par,sigma_par)) {
        return STLIM<T>::quiet_NaN();
    }

    //

    std::normal_distribution<T> norm_dist(T(0),T(1));

    return mu_par + sigma_par*norm_dist(engine);
}

template<typename T1, typename T2, typename TC = common_return_t<T1,T2>>
statslib_inline
TC
rnorm_type_check(const T1 mu_par, const T2 sigma_par, rand_engine_t& engine)
{
    return rnorm_compute(static_cast<TC>(mu_par),static_cast<TC>(sigma_par),engine);
}

}

/**
 * Random sampling function for the Normal distribution
 *
 * @param mu_par the mean parameter, a real-valued input.
 * @param sigma_par the standard deviation parameter, a real-valued input.
 * @param engine a random engine, passed by reference.
 *
 * @return a pseudo-random draw from the Normal distribution.
 *
 * Example:
 * \code{.cpp}
 * stats::rand_engine_t engine(1776);
 * stats::rnorm(1.0,2.0,engine);
 * \endcode
 */

template<typename T1, typename T2>
statslib_inline
common_return_t<T1,T2>
rnorm(const T1 mu_par, const T2 sigma_par, rand_engine_t& engine)
{
    return internal::rnorm_type_check(mu_par,sigma_par,engine);
}

/**
 * Random sampling function for the Normal distribution
 *
 * @param mu_par the mean parameter, a real-valued input.
 * @param sigma_par the standard deviation parameter, a real-valued input.
 * @param seed_val initialize the random engine with a non-negative integral-valued seed.
 *
 * @return a pseudo-random draw from the Normal distribution.
 *
 * Example:
 * \code{.cpp}
 * stats::rnorm(1.0,2.0,1776);
 * \endcode
 */

template<typename T1, typename T2>
statslib_inline
common_return_t<T1,T2>
rnorm(const T1 mu_par, const T2 sigma_par, const ullint_t seed_val)
{
    rand_engine_t engine(seed_val);
    return rnorm(mu_par,sigma_par,engine);
}

/**
 * Random sampling function for the Standard Normal distribution
 *
 * @return a pseudo-random draw from the Standard Normal distribution.
 *
 * Example:
 * \code{.cpp}
 * stats::rnorm();
 * \endcode
 */

template<typename T>
statslib_inline
T
rnorm()
{
    return rnorm(T(0),T(1));
}

//

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename T1, typename T2, typename rT>
statslib_inline
void
rnorm_vec(const T1 mu_par, const T2 sigma_par, rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    RAND_DIST_FN_VEC(rnorm,vals_out,num_elem,mu_par,sigma_par);
}
#endif

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename T2>
statslib_inline
void
rnorm_mat_check(std::vector<eT>& X, const T1 mu_par, const T2 sigma_par)
{
    STDVEC_RAND_DIST_FN(rnorm,mu_par,sigma_par);
}
#endif

#ifdef STATS_ENABLE_MATRIX_FEATURES
template<typename mT, typename T1, typename T2>
statslib_inline
void
rnorm_mat_check(mT& X, const T1 mu_par, const T2 sigma_par)
{
    MAIN_MAT_RAND_DIST_FN(rnorm,mu_par,sigma_par);
}
#endif

}

/**
 * @brief Random matrix sampling function for the Normal distribution
 *
 * @param n the number of output rows
 * @param k the number of output columns
 * @param mu_par the mean parameter, a real-valued input.
 * @param sigma_par the standard deviation parameter, a real-valued input.
 *
 * @return a matrix of pseudo-random draws from the Normal distribution.
 *
 * Example:
 * \code{.cpp}
 * // std::vector
 * stats::rnorm<std::vector<double>>(5,4,1.0,2.0);
 * // Armadillo matrix
 * stats::rnorm<arma::mat>(5,4,1.0,2.0);
 * // Blaze dynamic matrix
 * stats::rnorm<blaze::DynamicMatrix<double,blaze::columnMajor>>(5,4,1.0,2.0);
 * // Eigen dynamic matrix
 * stats::rnorm<Eigen::MatrixXd>(5,4,1.0,2.0);
 * \endcode
 *
 * @note This function requires template instantiation; acceptable output types include: <tt>std::vector</tt> with primitive types (e.g., \c float, \c double, etc.), as well as Armadillo, Blaze, and Eigen dense matrices.
 */

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename mT, typename T1, typename T2>
statslib_inline
mT
rnorm(const ullint_t n, const ullint_t k, const T1 mu_par, const T2 sigma_par)
{
    GEN_MAT_RAND_FN(rnorm_mat_check,mu_par,sigma_par);
}
#endif
