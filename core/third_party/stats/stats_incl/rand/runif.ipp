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
 * Sample from a uniform distribution
 */

//
// scalar output

namespace internal
{

template<typename T>
statslib_inline
T
runif_compute(const T a_par, const T b_par, rand_engine_t& engine)
{
    if (!unif_sanity_check(a_par,b_par)) {
        return STLIM<T>::quiet_NaN();
    }
    
    // convert from [a,b) to (a,b)

    T a_par_adj = std::nextafter(a_par, b_par);
    std::uniform_real_distribution<T> unif_dist(a_par_adj, b_par);

    return unif_dist(engine);
}

template<typename T1, typename T2, typename TC = common_return_t<T1,T2>>
statslib_inline
TC
runif_type_check(const T1 a_par, const T2 b_par, rand_engine_t& engine)
{
    return runif_compute(static_cast<TC>(a_par),static_cast<TC>(b_par),engine);
}

}

/**
 * @brief Random sampling function for the Uniform distribution
 *
 * @param a_par a real-valued shape parameter.
 * @param b_par a real-valued shape parameter.
 * @param engine a random engine, passed by reference.
 *
 * @return a pseudo-random draw from the Uniform distribution.
 *
 * Example:
 * \code{.cpp}
 * stats::rand_engine_t engine(1776);
 * stats::runif(3.0,2.0,engine);
 * \endcode
 */

template<typename T1, typename T2>
statslib_inline
common_return_t<T1,T2> 
runif(const T1 a_par, const T2 b_par, rand_engine_t& engine)
{
    return internal::runif_type_check(a_par,b_par,engine);
}

/**
 * @brief Random sampling function for the Uniform distribution
 *
 * @param a_par a real-valued shape parameter.
 * @param b_par a real-valued shape parameter.
 * @param seed_val initialize the random engine with a non-negative integral-valued seed.
 *
 * @return a pseudo-random draw from the Uniform distribution.
 *
 * Example:
 * \code{.cpp}
 * stats::runif(3.0,2.0,1776);
 * \endcode
 */

template<typename T1, typename T2>
statslib_inline
common_return_t<T1,T2> 
runif(const T1 a_par, const T2 b_par, const ullint_t seed_val)
{
    rand_engine_t engine(seed_val);
    return runif(a_par,b_par,engine);
}

/**
 * @brief Random sampling function for the Uniform distribution on the unit interval
 *
 * @return a pseudo-random draw from the Uniform distribution.
 *
 * Example:
 * \code{.cpp}
 * stats::runif();
 * \endcode
 */

template<typename T>
statslib_inline
T
runif()
{
    return runif(T(0),T(1));
}

//
// vector/matrix output

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename T1, typename T2, typename rT>
statslib_inline
void
runif_vec(const T1 a_par, const T2 b_par, rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    RAND_DIST_FN_VEC(runif,vals_out,num_elem,a_par,b_par);
}
#endif

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename T2>
statslib_inline
void
runif_mat_check(std::vector<eT>& X, const T1 a_par, const T2 b_par)
{
    STDVEC_RAND_DIST_FN(runif,a_par,b_par);
}
#endif

#ifdef STATS_ENABLE_MATRIX_FEATURES
template<typename mT, typename T1, typename T2>
statslib_inline
void
runif_mat_check(mT& X, const T1 a_par, const T2 b_par)
{
    MAIN_MAT_RAND_DIST_FN(runif,a_par,b_par);
}
#endif

}

/**
 * @brief Random matrix sampling function for the Uniform distribution
 *
 * @param n the number of output rows
 * @param k the number of output columns
 * @param a_par a real-valued shape parameter.
 * @param b_par a real-valued shape parameter.
 *
 * @return a matrix of pseudo-random draws from the Uniform distribution.
 *
 * Example:
 * \code{.cpp}
 * // std::vector
 * stats::runif<std::vector<double>>(5,4,-1.0,3.0);
 * // Armadillo matrix
 * stats::runif<arma::mat>(5,4,-1.0,3.0);
 * // Blaze dynamic matrix
 * stats::runif<blaze::DynamicMatrix<double,blaze::columnMajor>>(5,4,-1.0,3.0);
 * // Eigen dynamic matrix
 * stats::runif<Eigen::MatrixXd>(5,4,-1.0,3.0);
 * \endcode
 *
 * @note This function requires template instantiation; acceptable output types include: <tt>std::vector</tt> with primitive types (e.g., \c float, \c double, etc.), as well as Armadillo, Blaze, and Eigen dense matrices.
 */

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename mT, typename T1, typename T2>
statslib_inline
mT
runif(const ullint_t n, const ullint_t k, const T1 a_par, const T2 b_par)
{
    GEN_MAT_RAND_FN(runif_mat_check,a_par,b_par);
}
#endif
