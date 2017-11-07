/**  \file ublas_elementwise.hpp \brief Elementwise operations for ublas vector/matrix */
/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2015 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   BayesOpt is free software: you can redistribute it and/or modify it 
   under the terms of the GNU Affero General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   BayesOpt is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Affero General Public License for more details.

   You should have received a copy of the GNU Affero General Public License
   along with BayesOpt.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------
*/

#ifndef  _ELEMENTWISE_UBLAS_HPP_
#define  _ELEMENTWISE_UBLAS_HPP_

#include <algorithm>
#include <boost/numeric/ublas/vector.hpp>

namespace bayesopt
{
  namespace utils
  {

    /** 
     * Computes the elementwise product of two vectors or matrices.
     *             c_i = a_i * b_i
     */
    template <class v1, class v2>
    v1 ublas_elementwise_prod(const v1& a, const v2& b)
    {
      typedef typename v1::value_type D;
      v1 c(a.size());
      std::transform(a.begin(),a.end(),b.begin(),c.begin(),std::multiplies<D>());
      return c;
    }

    /** 
     * Computes the elementwise division of two vectors or matrices.
     *            c_i = a_i / b_i
     */
    template <class v1, class v2>
    v1 ublas_elementwise_div(const v1& a, const v2& b)
    {
      typedef typename v1::value_type D;
      v1 c(a.size());
      std::transform(a.begin(),a.end(),b.begin(),c.begin(),std::divides<D>());
      return c;
    }

  } //namespace utils
} //namespace bayesopt

#endif
