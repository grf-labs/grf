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

#ifndef  _KERNEL_RQ_HPP_
#define  _KERNEL_RQ_HPP_

#include "kernels/kernel_atomic.hpp"

namespace bayesopt
{
  
  /**\addtogroup KernelFunctions */
  //@{

  /** \brief Rational quadratic (Student's t) kernel. Isotropic version. */
  class RQIso: public ISOkernel
  {
  public:
    void init(size_t input_dim)
    { n_params = 2; n_inputs = input_dim; };

    double operator()( const vectord &x1, const vectord &x2)
    {
      double rl = computeWeightedNorm2(x1,x2);
      double k = rl*rl/(2*params(1));
      return std::pow(1+k,-params(1));
    };

    // TODO:
    double gradient(const vectord &x1, const vectord &x2,
		    size_t component)
    { assert(false); return 0.0;  };
  };

  //@}

} //namespace bayesopt

#endif
