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

#ifndef  _KERNEL_POLYNOMIAL_HPP_
#define  _KERNEL_POLYNOMIAL_HPP_

#include "kernels/kernel_atomic.hpp"

namespace bayesopt
{
  
  /**\addtogroup KernelFunctions */
  //@{

  /** Polynomial covariance function*/
  class Polynomial: public AtomicKernel
  {
  public:
    Polynomial(){ mExp = 1; };

    void init(size_t input_dim)
    { n_params = 2;  n_inputs = input_dim;  };

    double operator()( const vectord &x1, const vectord &x2)
    {
      double xx = boost::numeric::ublas::inner_prod(x1,x2); 
      return params(0)*params(0) * std::pow((params(1)+xx),static_cast<int>(mExp));
    };

    //TODO:
    double gradient( const vectord &x1, const vectord &x2,
		     size_t component)
    { assert(false); return 0.0; };
  protected:
    size_t mExp;
  };

  class Polynomial2: public Polynomial { public: Polynomial2(){ mExp = 2;};};
  class Polynomial3: public Polynomial { public: Polynomial3(){ mExp = 3;};};
  class Polynomial4: public Polynomial { public: Polynomial4(){ mExp = 4;};};
  class Polynomial5: public Polynomial { public: Polynomial5(){ mExp = 5;};};
  class Polynomial6: public Polynomial { public: Polynomial6(){ mExp = 6;};};
  class Polynomial7: public Polynomial { public: Polynomial7(){ mExp = 7;};};

  //@}

} //namespace bayesopt

#endif
