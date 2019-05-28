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
 * quantile function of the t distribution
 *
 * using   Hill, G. W. (1970) "Algorithm 396: Student's t-Quantiles", ACM COMS 13(10), 619-20
 *         Hill, G. W. (1981) "Remark on "Algorithm 396: Student's t-Quantiles"", ACM TOMS 7, 250-1
 */

//
// single input

// coefficients

namespace internal
{

template<typename T>
statslib_constexpr
T
qt_int_coef_a(const T dof_par)
noexcept
{
    return( T(1) / (dof_par - T(0.5)) );
}

template<typename T>
statslib_constexpr
T
qt_int_coef_b(const T coef_a)
noexcept
{
    return( T(48) / (coef_a * coef_a) );
}

template<typename T>
statslib_constexpr
T
qt_int_coef_c(const T coef_a, const T coef_b)
noexcept
{
    return( T(96.36) + coef_a * ( - T(16) + coef_a * ( - T(98) + T(20700)*coef_a/coef_b ) ) );
}

template<typename T>
statslib_constexpr
T
qt_int_coef_c_update(const T x, const T dof_par, const T coef_c)
noexcept
{
    return( dof_par < T(5) ? \
                coef_c + T(0.3)*(dof_par - T(4.5))*(x + T(0.6)) :
            // else
                coef_c );
}

template<typename T>
statslib_constexpr
T
qt_int_coef_d(const T dof_par, const T coef_a, const T coef_b, const T coef_c)
noexcept
{
    return( dof_par * stmath::sqrt(coef_a) * GCEM_SQRT_HALF_PI * \
                ( T(1) + (-T(3) + T(94.5)/(coef_b + coef_c))/coef_b ) );
}

// initial y

template<typename T>
statslib_constexpr
T
qt_int_y_init(const T p, const T dof_par, const T coef_d)
noexcept
{
    return( stmath::pow( coef_d*p , T(2)/dof_par ) );
}

// update y

template<typename T>
statslib_constexpr
T
qt_int_y_1_update(const int stage, const T y, const T x, const T coef_a, const T coef_b, const T coef_c)
noexcept
{
    return( stage == 1 ? \
                x*x :
            stage == 2 ? \
                x * ( T(1) + ( ( y * (T(36) + y * (T(6.3) + y*T(0.4))) + T(94.5) )/coef_c - y - T(3) ) / coef_b ) :
            stage == 3 ? \
                coef_a * y * y :
            //
            y > T(0.1) ? \
                stmath::expm1(y) : // see Remark on Algorithm 396
                y + y * y * (T(12) + y * (T(4) + y))/T(24) );
}

//

template<typename T>
statslib_constexpr
T
qt_int_y_1(const int stage, const T y, const T x, const T dof_par, const T coef_a, const T coef_b, const T coef_c, const T coef_d)
noexcept
{
    return( stage == 0 ? \
                qt_int_y_1(1,T(0), x, dof_par, coef_a, coef_b, qt_int_coef_c_update(x,dof_par,coef_c), coef_d) :
            //
            stage == 1 ? \
                qt_int_y_1(2, qt_int_y_1_update(1,y,x,coef_a,coef_b,coef_c), x, dof_par, coef_a, coef_b,
                           coef_b + coef_c + x * (- T(2) + x * (- T(7) + x * (-T(5) + T(0.05)*x*coef_d))), coef_d) :
            //
            stage == 2 ? \
                qt_int_y_1(3,qt_int_y_1_update(2,y,x,coef_a,coef_b,coef_c),x,dof_par,coef_a,coef_b,coef_c,coef_d) :
            stage == 3 ? \
                qt_int_y_1(4,qt_int_y_1_update(3,y,x,coef_a,coef_b,coef_c),x,dof_par,coef_a,coef_b,coef_c,coef_d) :
            //
                qt_int_y_1_update(4,y,x,coef_a,coef_b,coef_c) );
}

// small dof cases

template<typename T>
statslib_constexpr
T qt_int_y_2_iter_loop(const T p, const T dof_par, const T val_out, const T x, const int iter) noexcept;

template<typename T>
statslib_constexpr
T
qt_int_y_2_iter_val_update(const T val_out, const T dof_par, const T x)
noexcept
{
    return val_out + x * (T(1) + x*val_out * (T(1) + dof_par) / (T(2) * (dof_par + val_out*val_out)));
}

template<typename T>
statslib_constexpr
T
qt_int_y_2_iter_loop_update(const T p, const T dof_par, const T val_out, const int iter)
noexcept
{
    return( iter < 6 ? \
                qt_int_y_2_iter_loop(p,dof_par,val_out,(T(1)-pt(val_out,dof_par,false)-p)/dt(val_out,dof_par,false),iter+1) :
            // else
                val_out );
}

template<typename T>
statslib_constexpr
T
qt_int_y_2_iter_loop(const T p, const T dof_par, const T val_out, const T x, const int iter)
noexcept
{
    return qt_int_y_2_iter_loop_update(p,dof_par,qt_int_y_2_iter_val_update(val_out,dof_par,x),iter);
}

template<typename T>
statslib_constexpr
T
qt_int_y_2_iter_begin(const T p, const T dof_par, const T y)
noexcept
{
    return qt_int_y_2_iter_loop(p/T(2),dof_par,stmath::sqrt(dof_par*y),T(0),0);
}

template<typename T>
statslib_constexpr
T
qt_int_y_2_mterm(const T y, const T dof_par, const T coef_d)
noexcept
{
    return( T(0.5) / (dof_par + T(4)) + T(1) / ( ((dof_par + T(6)) / (dof_par * y) - T(0.089)*coef_d - T(0.822)) * (dof_par + T(2)) * T(3) ) );
}

template<typename T>
statslib_constexpr
T
qt_int_y_2(const T p, const T y, const T dof_par, const T coef_d)
noexcept
{
    return qt_int_y_2_iter_begin(p, dof_par, (T(1)/y) + (dof_par + T(1)) * ( qt_int_y_2_mterm(y,dof_par,coef_d) * y - T(1)) / (dof_par + T(2)));
}

//

template<typename T>
statslib_constexpr
T
qt_int_choose(const T p, const T y, const T dof_par, const T coef_a, const T coef_b, const T coef_c, const T coef_d)
noexcept
{
    return( y > T(0.05) + coef_a ? \
                stmath::sqrt( dof_par * qt_int_y_1(0,y,qnorm(T(0.5)*p),dof_par,coef_a,coef_b,coef_c,coef_d)) :
            // else
                qt_int_y_2(p,y,dof_par,coef_d) );
}

template<typename T>
statslib_constexpr
T
qt_int_finish(const T p, const T dof_par, const T coef_a, const T coef_b, const T coef_c, const T coef_d)
noexcept
{
    return qt_int_choose(p, qt_int_y_init(p,dof_par,coef_d), dof_par,coef_a,coef_b,coef_c,coef_d);
}

template<typename T>
statslib_constexpr
T
qt_int_main_iter(const ullint_t stage, const T p, const T dof_par, const T coef_a, const T coef_b, const T coef_c, const T coef_d)
noexcept
{
    return( stage == 0U ? \
                qt_int_main_iter(1U,p,dof_par,qt_int_coef_a(dof_par),T(0),T(0),T(0)) :
            stage == 1U ? \
                qt_int_main_iter(2U,p,dof_par,coef_a,qt_int_coef_b(coef_a),T(0),T(0)) :
            stage == 2U ? \
                qt_int_main_iter(3U,p,dof_par,coef_a,coef_b,qt_int_coef_c(coef_a,coef_b),T(0)) :
            stage == 3U ? \
                qt_int_main_iter(4U,p,dof_par,coef_a,coef_b,coef_c,qt_int_coef_d(dof_par,coef_a,coef_b,coef_c)) :
            // else
                qt_int_finish(p,dof_par,coef_a,coef_b,coef_c,coef_d) );
}

template<typename T>
statslib_constexpr
T
qt_int_main(const T p, const T dof_par)
noexcept
{
    return( p < T(0.5) ? \
                - qt_int_main_iter(0U,2*p,dof_par,T(0),T(0),T(0),T(0)) : 
                  qt_int_main_iter(0U,2*p,dof_par,T(0),T(0),T(0),T(0)) );
}

template<typename T>
statslib_constexpr
T
qt_vals_check(const T p, const T dof_par)
noexcept
{
    return( !t_sanity_check(dof_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            !prob_val_check(p) ? \
                STLIM<T>::quiet_NaN() :
            //
            p == T(0) ? \
                - STLIM<T>::infinity() :
            p == T(1) ? \
                STLIM<T>::infinity() :
            // Cauchy case, etc.
            dof_par == T(1) ? \
                stmath::tan(GCEM_PI*(p - T(0.5))) :
            dof_par == T(2) ? \
                (2*p - T(1)) / stmath::sqrt(2*p*(T(1) - p)) :
            // normal case
            dof_par == STLIM<T>::infinity() ? \
                qnorm(p,T(0),T(1)) :
            // else
                qt_int_main(p,dof_par) );
}

template<typename T1, typename T2, typename TC = common_return_t<T1,T2>>
statslib_constexpr
TC
qt_type_check(const T1 p, const T2 dof_par)
noexcept
{
    return qt_vals_check(static_cast<TC>(p),static_cast<TC>(dof_par));
}

}

/**
 * @brief Quantile function of the t-distribution
 *
 * @param p a real-valued input.
 * @param dof_par the degrees of freedom parameter, a real-valued input.
 *
 * @return the quantile function evaluated at \c p.
 * 
 * Example:
 * \code{.cpp} stats::qt(0.5,11); \endcode
 */

template<typename T1, typename T2>
statslib_constexpr
common_return_t<T1,T2>
qt(const T1 p, const T2 dof_par)
noexcept
{
    return internal::qt_type_check(p,dof_par);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename rT>
statslib_inline
void
qt_vec(const eT* __stats_pointer_settings__ vals_in, const T1 dof_par, 
             rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(qt,vals_in,vals_out,num_elem,dof_par);
}
#endif

}

/**
 * @brief Quantile function of the t-distribution
 *
 * @param x a standard vector.
 * @param dof_par the degrees of freedom parameter, a real-valued input.
 *
 * @return a vector of quantile values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<double> x = {0.3, 0.5, 0.8};
 * stats::qt(x,4);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
std::vector<rT>
qt(const std::vector<eT>& x, const T1 dof_par)
{
    STDVEC_DIST_FN(qt_vec,dof_par);
}
#endif

/**
 * @brief Quantile function of the t-distribution
 *
 * @param X a matrix of input values.
 * @param dof_par the degrees of freedom parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * arma::mat X = { {0.2, 0.7, 0.9},
 *                 {0.1, 0.8, 0.3} };
 * stats::qt(X,4);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
ArmaMat<rT>
qt(const ArmaMat<eT>& X, const T1 dof_par)
{
    ARMA_DIST_FN(qt_vec,dof_par);
}

template<typename mT, typename tT, typename T1>
statslib_inline
mT
qt(const ArmaGen<mT,tT>& X, const T1 dof_par)
{
    return qt(X.eval(),dof_par);
}
#endif

/**
 * @brief Quantile function of the t-distribution
 *
 * @param X a matrix of input values.
 * @param dof_par the degrees of freedom parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::qt(X,4);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
qt(const BlazeMat<eT,To>& X, const T1 dof_par)
{
    BLAZE_DIST_FN(qt_vec,dof_par);
}
#endif

/**
 * @brief Quantile function of the t-distribution
 *
 * @param X a matrix of input values.
 * @param dof_par the degrees of freedom parameter, a real-valued input.
 *
 * @return a matrix of quantile values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::qt(X,4);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
qt(const EigenMat<eT,iTr,iTc>& X, const T1 dof_par)
{
    EIGEN_DIST_FN(qt_vec,dof_par);
}
#endif
