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
 * Sample from a Weibull distribution
 */

//
// scalar output

namespace internal
{

template<typename T>
statslib_inline
T
rweibull_compute(const T shape_par, const T scale_par, rand_engine_t& engine)
{
    return( !weibull_sanity_check(shape_par,scale_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            qweibull(runif(T(0),T(1),engine),shape_par,scale_par) );
}

template<typename T1, typename T2, typename TC = common_return_t<T1,T2>>
statslib_constexpr
TC
rweibull_type_check(const T1 shape_par, const T2 scale_par, rand_engine_t& engine)
noexcept
{
    return rweibull_compute(static_cast<TC>(shape_par),static_cast<TC>(scale_par),engine);
}

}

/**
 * @brief Random sampling function for the Weibull distribution
 *
 * @param shape_par the shape parameter, a real-valued input.
 * @param scale_par the scale parameter, a real-valued input.
 * @param engine a random engine, passed by reference.
 *
 * @return a pseudo-random draw from the Weibull distribution.
 *
 * Example:
 * \code{.cpp}
 * stats::rand_engine_t engine(1776);
 * stats::rweibull(3.0,2.0,engine);
 * \endcode
 */

template<typename T1, typename T2>
statslib_inline
common_return_t<T1,T2>
rweibull(const T1 shape_par, const T2 scale_par, rand_engine_t& engine)
{
    return internal::rweibull_type_check(shape_par,scale_par,engine);
}

/**
 * @brief Random sampling function for the Weibull distribution
 *
 * @param shape_par the shape parameter, a real-valued input.
 * @param scale_par the scale parameter, a real-valued input.
 * @param seed_val initialize the random engine with a non-negative integral-valued seed.
 *
 * @return a pseudo-random draw from the Weibull distribution.
 *
 * Example:
 * \code{.cpp}
 * stats::rweibull(3.0,2.0,1776);
 * \endcode
 */

template<typename T1, typename T2>
statslib_inline
common_return_t<T1,T2>
rweibull(const T1 shape_par, const T2 scale_par, const ullint_t seed_val)
{
    rand_engine_t engine(seed_val);
    return rweibull(shape_par,scale_par,engine);
}

//
// vector/matrix output

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename T1, typename T2, typename rT>
statslib_inline
void
rweibull_vec(const T1 shape_par, const T2 scale_par, rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    RAND_DIST_FN_VEC(rweibull,vals_out,num_elem,shape_par,scale_par);
}
#endif

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename T2>
statslib_inline
void
rweibull_mat_check(std::vector<eT>& X, const T1 shape_par, const T2 scale_par)
{
    STDVEC_RAND_DIST_FN(rweibull,shape_par,scale_par);
}
#endif

#ifdef STATS_ENABLE_MATRIX_FEATURES
template<typename mT, typename T1, typename T2>
statslib_inline
void
rweibull_mat_check(mT& X, const T1 shape_par, const T2 scale_par)
{
    MAIN_MAT_RAND_DIST_FN(rweibull,shape_par,scale_par);
}
#endif

}

/**
 * @brief Random matrix sampling function for the Weibull distribution
 *
 * @param n the number of output rows
 * @param k the number of output columns
 * @param shape_par the shape parameter, a real-valued input.
 * @param scale_par the scale parameter, a real-valued input.
 *
 * @return a matrix of pseudo-random draws from the Weibull distribution.
 *
 * Example:
 * \code{.cpp}
 * // std::vector
 * stats::rweibull<std::vector<double>>(5,4,3.0,2.0);
 * // Armadillo matrix
 * stats::rweibull<arma::mat>(5,4,3.0,2.0);
 * // Blaze dynamic matrix
 * stats::rweibull<blaze::DynamicMatrix<double,blaze::columnMajor>>(5,4,3.0,2.0);
 * // Eigen dynamic matrix
 * stats::rweibull<Eigen::MatrixXd>(5,4,3.0,2.0);
 * \endcode
 *
 * @note This function requires template instantiation; acceptable output types include: <tt>std::vector</tt> with primitive types (e.g., \c float, \c double, etc.), as well as Armadillo, Blaze, and Eigen dense matrices.
 */

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename mT, typename T1, typename T2>
statslib_inline
mT
rweibull(const ullint_t n, const ullint_t k, const T1 shape_par, const T2 scale_par)
{
    GEN_MAT_RAND_FN(rweibull_mat_check,shape_par,scale_par);
}
#endif
