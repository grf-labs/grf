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

#ifndef  _KERNEL_HAMMING_HPP_
#define  _KERNEL_HAMMING_HPP_

#include "kernels/kernel_atomic.hpp"

namespace bayesopt
{
  
  /**\addtogroup KernelFunctions */
  //@{

  /** Kernel for categorical data. It measures the hamming distance
   *  between vectors. 
   */
  class HammingKernel: public AtomicKernel
  {
  public:
    void init(size_t input_dim)
    { n_params = 1; n_inputs = input_dim;  };

    size_t hammingDistance(const vectori& s1, const vectori& s2)
    {
      size_t hdist = 0;
      vectori::const_iterator i1;
      vectori::const_iterator i2;
      for( i1 = s1.begin(), i2 = s2.begin(); 
	   i1 < s1.end() && i2 < s2.end();
	   ++i1, ++i2 )
	{
	  hdist += (*i1 == *i2) ? 0 : 1;
	}
      return hdist;
    }

    double operator()(const vectord &x1, const vectord &x2)
    { 
      const size_t n = x1.size();
      const double coef = -params(0)/2.0;
      vectori s1(n);
      vectori s2(n);

      for(size_t i=0; i<n; ++i)
	{
	  // We add 0.5 to avoid floating point approximation errors
	  s1(i) = static_cast<int>(x1(i)+0.5);
	  s2(i) = static_cast<int>(x2(i)+0.5);
	}
      
      const double dist = static_cast<double>(hammingDistance(s1,s2));
      return std::exp(coef*dist*dist);
    };

    double gradient(const vectord &x1, const vectord &x2,
		    size_t component)
    { return 0.0; };
  };

  //@}

} //namespace bayesopt

#endif
