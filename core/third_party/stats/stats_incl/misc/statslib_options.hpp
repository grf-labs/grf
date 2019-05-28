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

#ifndef _statslib_options_HPP
#define _statslib_options_HPP

//
// StatsLib options and compiler directives

// version

#ifndef STATS_VERSION_MAJOR
    #define STATS_VERSION_MAJOR 3
#endif

#ifndef STATS_VERSION_MINOR
    #define STATS_VERSION_MINOR 0
#endif

#ifndef STATS_VERSION_PATCH
    #define STATS_VERSION_PATCH 0
#endif

// switch between inline mode vs constexpr

#ifndef statslib_inline
    #define statslib_inline inline
#endif

#ifndef STATS_GO_INLINE
    #define statslib_constexpr constexpr
    #define stmath gcem
#else
    #define statslib_constexpr inline
    #include <cmath>
    #define stmath std
#endif

// include some basic libraries

#include <limits>
#include <random>

// typedefs

namespace stats
{
    using uint_t = unsigned int;
    using ullint_t = unsigned long long int;

    using llint_t = long long int;

    using rand_engine_t = std::mt19937_64;

    namespace GCINT = gcem::internal;

    template<class T>
    using STLIM = std::numeric_limits<T>;

    template<typename T>
    using return_t = typename std::conditional<std::is_integral<T>::value,double,T>::type;

    template<typename ...T>
    using common_t = typename std::common_type<T...>::type;

    template<typename ...T>
    using common_return_t = return_t<common_t<T...>>;
}

// enable OpenMP

#if defined(_OPENMP) && !defined(STATS_DONT_USE_OPENMP) && !defined(STATS_USE_OPENMP)
    #define STATS_USE_OPENMP
    #include <omp.h>

    #define STATS_OMP_N_BLOCKS ullint_t(4)
    #define STATS_OMP_MIN_N_PER_BLOCK ullint_t(4)
#endif

// enable std::vector features

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
    #include <algorithm>
    #include <vector>
#endif

// #ifndef STATS_ENABLE_STDVEC_WRAPPERS
//     #if defined(_LIBCPP_VECTOR) || defined(_GLIBCXX_VECTOR)
//         #define STATS_ENABLE_STDVEC_WRAPPERS
//     #endif
// #endif

// enable wrappers for linear algebra libraries

#if defined(STATS_ENABLE_ARMA_WRAPPERS) || defined(STATS_ENABLE_BLAZE_WRAPPERS) || defined(STATS_ENABLE_EIGEN_WRAPPERS)
    #define STATS_ENABLE_MATRIX_FEATURES
#endif

//

#if defined(STATS_ENABLE_STDVEC_WRAPPERS) || defined(STATS_ENABLE_MATRIX_FEATURES)
    #define STATS_ENABLE_INTERNAL_VEC_FEATURES
#endif

// Armadillo options

#ifdef STATS_ENABLE_ARMA_WRAPPERS
    #ifdef USE_RCPP_ARMADILLO
        #include <RcppArmadillo.h>
    #else
        #ifndef ARMA_DONT_USE_WRAPPER
            #define ARMA_DONT_USE_WRAPPER
        #endif
        #include "armadillo"
    #endif

    template<typename eT>
    using ArmaMat = arma::Mat<eT>;

    template<typename mT, typename tT>
    using ArmaGen = arma::Gen<mT,tT>;

    template<typename T>
    using not_arma_mat = std::enable_if<!(std::is_same<T,arma::Mat<double>>::value || \
                                          std::is_same<T,arma::Mat<float>>::value)>;
    
    #if defined(STATS_ENABLE_BLAZE_WRAPPERS) || defined(STATS_ENABLE_EIGEN_WRAPPERS)
        #error StatsLib cannot interface with more than one matrix library at a time
    #endif
#else
    template<typename T>
    using not_arma_mat = std::enable_if<true>;
#endif

// Blaze options

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
    #include "blaze/Blaze.h"
    #include <iostream>

    template<typename eT, bool To = blaze::columnMajor>
    using BlazeMat = blaze::DynamicMatrix<eT,To>;

    template<typename eT, bool To = blaze::rowMajor>
    using BlazeRow = blaze::DynamicVector<eT,To>;

    template<typename T>
    using not_blaze_mat = std::enable_if<!(std::is_same<T,BlazeMat<double,blaze::columnMajor>>::value || \
                                           std::is_same<T,BlazeMat<double,blaze::rowMajor>>::value ||    \
                                           std::is_same<T,BlazeMat<float,blaze::columnMajor>>::value ||  \
                                           std::is_same<T,BlazeMat<float,blaze::rowMajor>>::value)>;
#endif

// Eigen Options

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
    #include <Eigen/Dense>
    #include <iostream>

    template<typename eT, int iTr, int iTc>
    using EigenMat = Eigen::Matrix<eT,iTr,iTc>;

    #if defined(STATS_ENABLE_BLAZE_WRAPPERS)
        #error StatsLib cannot interface with more than one matrix library at a time
    #endif

#endif

// other

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
    #include <iostream> // cerr, cout, ...
#endif

//
// misc. compiler options

#ifndef __stats_pointer_settings__
    #if defined(__clang__) || defined(__GNUC__)
        #define __stats_pointer_settings__ __restrict__
    #elif defined(_MSC_VER)
        #define __stats_pointer_settings__ __restrict
    #else
        #define __stats_pointer_settings__
    #endif
#endif

#endif
