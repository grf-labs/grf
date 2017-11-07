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

#ifndef  _KERNEL_PROD_HPP_
#define  _KERNEL_PROD_HPP_

#include "kernels/kernel_combined.hpp"

namespace bayesopt
{
  
  /**\addtogroup KernelFunctions */
  //@{


  /** \brief Product of two kernels */
  class KernelProd: public CombinedKernel
  {
  public:
    double operator()(const vectord &x1, const vectord &x2)
    { return (*left)(x1,x2) * (*right)(x1,x2); };

    //TODO: Not implemented
    double gradient(const vectord &x1, const vectord &x2,
		    size_t component)
    { return 0.0; };
  };

  //@}

} //namespace bayesopt

#endif
