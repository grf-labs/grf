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
 * Sample from a gamma distribution
 * using the Marsaglia and Tsang method
 */

//
// scalar output

namespace internal
{

template<typename T>
statslib_inline
T
rgamma_compute(const T shape_par, const T scale_par, rand_engine_t& engine)
{
    if (!gamma_sanity_check(shape_par,scale_par)) {
        return STLIM<T>::quiet_NaN();
    }

    //

    T ret = 0;

    if (shape_par > T(1))
    {
        const T d = shape_par - T(1)/T(3);
        const T c = (T(1) / T(3)) / std::sqrt(d);
        T V = 1;

        bool keep_running = true;

        while (keep_running)
        {
            T Z = rnorm(T(0),T(1),engine);

            if (Z > -T(1)/c)
            {
                V = std::pow(T(1) + c*Z,3);
                T U = runif(T(0),T(1),engine);

                T check_2 = T(0.5)*Z*Z + d*(T(1) - V + std::log(V));

                if (std::log(U) < check_2) {
                    keep_running = false;
                }
            }
        }

        ret = d * V * scale_par;
    }
    else
    {
        const T U = runif(T(0),T(1),engine);
        ret = rgamma(T(1) + shape_par, scale_par,engine) * std::pow(U,T(1)/shape_par);
    }

    //

    return ret;
}

template<typename T1, typename T2, typename TC = common_return_t<T1,T2>>
statslib_constexpr
TC
rgamma_type_check(const T1 shape_par, const T2 scale_par, rand_engine_t& engine)
noexcept
{
    return rgamma_compute(static_cast<TC>(shape_par),static_cast<TC>(scale_par),engine);
}

}

/**
 * @brief Random sampling function for the Gamma distribution
 *
 * @param shape_par the shape parameter, a real-valued input.
 * @param scale_par the scale parameter, a real-valued input.
 * @param engine a random engine, passed by reference.
 *
 * @return a pseudo-random draw from the Gamma distribution.
 *
 * Example:
 * \code{.cpp}
 * stats::rand_engine_t engine(1776);
 * stats::rgamma(3.0,2.0,engine);
 * \endcode
 */

template<typename T1, typename T2>
statslib_inline
common_return_t<T1,T2>
rgamma(const T1 shape_par, const T2 scale_par, rand_engine_t& engine)
{
    return internal::rgamma_type_check(shape_par,scale_par,engine);
}

/**
 * @brief Random sampling function for the Gamma distribution
 *
 * @param shape_par the shape parameter, a real-valued input.
 * @param scale_par the scale parameter, a real-valued input.
 * @param seed_val initialize the random engine with a non-negative integral-valued seed.
 *
 * @return a pseudo-random draw from the Gamma distribution.
 *
 * Example:
 * \code{.cpp}
 * stats::rgamma(3.0,2.0,1776);
 * \endcode
 */

template<typename T1, typename T2>
statslib_inline
common_return_t<T1,T2>
rgamma(const T1 shape_par, const T2 scale_par, ullint_t seed_val)
{
    rand_engine_t engine(seed_val);
    return rgamma(shape_par,scale_par,engine);
}

//
// vector/matrix output

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename T1, typename T2, typename rT>
statslib_inline
void
rgamma_vec(const T1 shape_par, const T2 scale_par, rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    RAND_DIST_FN_VEC(rgamma,vals_out,num_elem,shape_par,scale_par);
}
#endif

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename T2>
statslib_inline
void
rgamma_mat_check(std::vector<eT>& X, const T1 shape_par, const T2 scale_par)
{
    STDVEC_RAND_DIST_FN(rgamma,shape_par,scale_par);
}
#endif

#ifdef STATS_ENABLE_MATRIX_FEATURES
template<typename mT, typename T1, typename T2>
statslib_inline
void
rgamma_mat_check(mT& X, const T1 shape_par, const T2 scale_par)
{
    MAIN_MAT_RAND_DIST_FN(rgamma,shape_par,scale_par);
}
#endif

}

/**
 * @brief Random matrix sampling function for the Gamma distribution
 *
 * @param n the number of output rows
 * @param k the number of output columns
 * @param shape_par the shape parameter, a real-valued input.
 * @param scale_par the scale parameter, a real-valued input.
 *
 * @return a matrix of pseudo-random draws from the Gamma distribution.
 *
 * Example:
 * \code{.cpp}
 * // std::vector
 * stats::rgamma<std::vector<double>>(5,4,3.0,2.0);
 * // Armadillo matrix
 * stats::rgamma<arma::mat>(5,4,3.0,2.0);
 * // Blaze dynamic matrix
 * stats::rgamma<blaze::DynamicMatrix<double,blaze::columnMajor>>(5,4,3.0,2.0);
 * // Eigen dynamic matrix
 * stats::rgamma<Eigen::MatrixXd>(5,4,3.0,2.0);
 * \endcode
 *
 * @note This function requires template instantiation; acceptable output types include: <tt>std::vector</tt> with primitive types (e.g., \c float, \c double, etc.), as well as Armadillo, Blaze, and Eigen dense matrices.
 */

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename mT, typename T1, typename T2>
statslib_inline
mT
rgamma(const ullint_t n, const ullint_t k, const T1 shape_par, const T2 scale_par)
{
    GEN_MAT_RAND_FN(rgamma_mat_check,shape_par,scale_par);
}
#endif
