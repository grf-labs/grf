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

#ifndef  _KERNEL_LINEAR_HPP_
#define  _KERNEL_LINEAR_HPP_

#include "kernels/kernel_atomic.hpp"

namespace bayesopt
{
  
  /**\addtogroup KernelFunctions */
  //@{

  /** \brief Linear kernel. */
  class LinKernel: public AtomicKernel
  {
  public:
    void init(size_t input_dim)
    { n_params = 0;  n_inputs = input_dim; };

    double operator()(const vectord &x1, const vectord &x2)
    {
      assert(x1.size() == x2.size());
      return boost::numeric::ublas::inner_prod(x1,x2); 
    };

    // TODO:
    double gradient(const vectord &x1, const vectord &x2,
		    size_t component)
    { assert(false);  return 0.0;   };
  };

  /** \brief Linear kernel. */
  class LinKernelARD: public AtomicKernel
  {
  public:
    void init(size_t input_dim)
    { n_params = input_dim;  n_inputs = input_dim;  };

    double operator()(const vectord &x1, const vectord &x2)
    {
      assert(x1.size() == x2.size());
      vectord v1 = utils::ublas_elementwise_div(x1, params);
      vectord v2 = utils::ublas_elementwise_div(x2, params);
      return boost::numeric::ublas::inner_prod(v1,v2); 
    };

    // TODO:
    double gradient(const vectord &x1, const vectord &x2,
		    size_t component)
    { assert(false); return 0.0;  };
  };

  //@}

} //namespace bayesopt

#endif
