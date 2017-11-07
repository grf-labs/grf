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

#ifndef  _KERNEL_MATERN_HPP_
#define  _KERNEL_MATERN_HPP_

#include "kernels/kernel_atomic.hpp"

namespace bayesopt
{
  
  /**\addtogroup KernelFunctions */
  //@{


  /** \brief Matern isotropic kernel of 1st order */
  class MaternIso1: public ISOkernel
  {
  public:
    void init(size_t input_dim)
    { n_params = 1;  n_inputs = input_dim;  };

    double operator()(const vectord &x1, const vectord &x2)
    {
      double r = computeWeightedNorm2(x1,x2);
      return exp(-r);
    };

    double gradient(const vectord &x1, const vectord &x2,
		    size_t component)
    {
      double r = computeWeightedNorm2(x1,x2);
      return r*exp(-r);
    };
  };


  /** \brief Matern ARD kernel of 1st order */
  class MaternARD1: public ARDkernel
  {
  public:
    void init(size_t input_dim)
    { n_params = input_dim; n_inputs = input_dim;  };

    double operator()(const vectord &x1, const vectord &x2)
    {
      double r = computeWeightedNorm2(x1,x2);
      return exp(-r);
    };

    //TODO: 
    double gradient(const vectord &x1, const vectord &x2,
		    size_t component)
    { assert(false);  return 0.0;  };
  };


  /** \brief Matern kernel of 3rd order */
  class MaternIso3: public ISOkernel
  {
  public:
    void init(size_t input_dim)
    { n_params = 1; n_inputs = input_dim;  };

    double operator()( const vectord &x1, const vectord &x2)
    {
      double r = sqrt(3.0) * computeWeightedNorm2(x1,x2);
      double er = exp(-r);
      return (1+r)*er;
    };

    double gradient( const vectord &x1, const vectord &x2,
		     size_t component)
    {
      double r = sqrt(3.0) * computeWeightedNorm2(x1,x2);
      double er = exp(-r);
      return r*r*er; 
    };
  };

  /** \brief Matern ARD kernel of 3rd order */
  class MaternARD3: public ARDkernel
  {
  public:
    void init(size_t input_dim)
    { n_params = input_dim;  n_inputs = input_dim; };

    double operator()( const vectord &x1, const vectord &x2)
    {
      double r = sqrt(3.0) * computeWeightedNorm2(x1,x2);
      double er = exp(-r);
      return (1+r)*er;
    };

    double gradient( const vectord &x1, const vectord &x2,
		     size_t component)
    {
      assert(false); return 0.0;
    };
  };


  /** \brief Matern isotropic kernel of 5th order */
  class MaternIso5: public ISOkernel
  {
  public:
    void init(size_t input_dim)
    { n_params = 1; n_inputs = input_dim;  };

    double operator()( const vectord &x1, const vectord &x2)
    {
      double r = sqrt(5.0) * computeWeightedNorm2(x1,x2);
      double er = exp(-r);
      return (1+r*(1+r/3))*er;
    };
    double gradient( const vectord &x1, const vectord &x2,
		     size_t component)
    {    
      double r = sqrt(5.0) * computeWeightedNorm2(x1,x2);
      double er = exp(-r);
      return r*(1+r)/3*r*er; 
    };
  };


  /** \brief Matern ARD kernel of 5th order */
  class MaternARD5: public ARDkernel
  {
  public:
    void init(size_t input_dim)
    { n_params = input_dim;  n_inputs = input_dim; };

    double operator()( const vectord &x1, const vectord &x2)
    {
      double r = sqrt(5.0) * computeWeightedNorm2(x1,x2);
      double er = exp(-r);
      return (1+r*(1+r/3))*er;
    };

    //TODO:
    double gradient( const vectord &x1, const vectord &x2,
		     size_t component)
    { assert(false); return 0.0; };
  };

  //@}

} //namespace bayesopt

#endif
