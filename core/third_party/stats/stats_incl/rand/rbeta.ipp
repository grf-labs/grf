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
 * Sample from a Beta distribution
 */

//
// scalar output

namespace internal
{

template<typename T>
statslib_inline
T
rbeta_compute(const T a_par, const T b_par, rand_engine_t& engine)
{
    if (!beta_sanity_check(a_par,b_par)) {
        return STLIM<T>::quiet_NaN();
    }

    //
    
    const T X = rgamma(a_par,T(1),engine);
    const T Y = rgamma(b_par,T(1),engine);
    
    return X / (X + Y);
}

template<typename T1, typename T2, typename TC = common_return_t<T1,T2>>
statslib_inline
TC
rbeta_type_check(const T1 a_par, const T2 b_par, rand_engine_t& engine)
{
    return rbeta_compute(static_cast<TC>(a_par),static_cast<TC>(b_par),engine);
}

}

/**
 * @brief Random sampling function for the Beta distribution
 *
 * @param a_par a real-valued shape parameter.
 * @param b_par a real-valued shape parameter.
 * @param engine a random engine, passed by reference.
 *
 * @return a pseudo-random draw from the Beta distribution.
 *
 * Example:
 * \code{.cpp}
 * stats::rand_engine_t engine(1776);
 * stats::rbeta(3.0,2.0,engine);
 * \endcode
 */

template<typename T1, typename T2>
statslib_inline
common_return_t<T1,T2> 
rbeta(const T1 a_par, const T2 b_par, rand_engine_t& engine)
{
    return internal::rbeta_type_check(a_par,b_par,engine);
}

/**
 * @brief Random sampling function for the Beta distribution
 *
 * @param a_par a real-valued shape parameter.
 * @param b_par a real-valued shape parameter.
 * @param seed_val initialize the random engine with a non-negative integral-valued seed.
 *
 * @return a pseudo-random draw from the Beta distribution.
 *
 * Example:
 * \code{.cpp}
 * stats::rbeta(3.0,2.0,1776);
 * \endcode
 */

template<typename T1, typename T2>
statslib_inline
common_return_t<T1,T2> 
rbeta(const T1 a_par, const T2 b_par, const ullint_t seed_val)
{
    rand_engine_t engine(seed_val);
    return rbeta(a_par,b_par,engine);
}

//
// vector/matrix output

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename T1, typename T2, typename rT>
statslib_inline
void
rbeta_vec(const T1 a_par, const T2 b_par, rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    RAND_DIST_FN_VEC(rbeta,vals_out,num_elem,a_par,b_par);
}
#endif

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename T2>
statslib_inline
void
rbeta_mat_check(std::vector<eT>& X, const T1 a_par, const T2 b_par)
{
    STDVEC_RAND_DIST_FN(rbeta,a_par,b_par);
}
#endif

#ifdef STATS_ENABLE_MATRIX_FEATURES
template<typename mT, typename T1, typename T2>
statslib_inline
void
rbeta_mat_check(mT& X, const T1 a_par, const T2 b_par)
{
    MAIN_MAT_RAND_DIST_FN(rbeta,a_par,b_par);
}
#endif

}

/**
 * @brief Random matrix sampling function for the Beta distribution
 *
 * @param n the number of output rows
 * @param k the number of output columns
 * @param a_par a real-valued shape parameter.
 * @param b_par a real-valued shape parameter.
 *
 * @return a matrix of pseudo-random draws from the Beta distribution.
 *
 * Example:
 * \code{.cpp}
 * // std::vector
 * stats::rbeta<std::vector<double>>(5,4,3.0,2.0);
 * // Armadillo matrix
 * stats::rbeta<arma::mat>(5,4,3.0,2.0);
 * // Blaze dynamic matrix
 * stats::rbeta<blaze::DynamicMatrix<double,blaze::columnMajor>>(5,4,3.0,2.0);
 * // Eigen dynamic matrix
 * stats::rbeta<Eigen::MatrixXd>(5,4,3.0,2.0);
 * \endcode
 *
 * @note This function requires template instantiation; acceptable output types include: <tt>std::vector</tt> with primitive types (e.g., \c float, \c double, etc.), as well as Armadillo, Blaze, and Eigen dense matrices.
 */

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename mT, typename T1, typename T2>
statslib_inline
mT
rbeta(const ullint_t n, const ullint_t k, const T1 a_par, const T2 b_par)
{
    GEN_MAT_RAND_FN(rbeta_mat_check,a_par,b_par);
}
#endif
