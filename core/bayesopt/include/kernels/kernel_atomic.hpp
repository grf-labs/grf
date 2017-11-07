/** \file kernel_atomic.hpp \brief Atomic (simple) kernel functions */
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

#ifndef  _KERNEL_ATOMIC_HPP_
#define  _KERNEL_ATOMIC_HPP_

#include <valarray>
#include "kernel_functors.hpp"
#include "ublas_elementwise.hpp"

namespace bayesopt
{
  
  /**\addtogroup KernelFunctions */
  //@{


  /** \brief Abstract class for an atomic kernel */
  class AtomicKernel : public Kernel
  {
  public:
    virtual void init(size_t input_dim)
    { n_inputs = input_dim; };

    void setHyperParameters(const vectord &theta) 
    {
      if(theta.size() != n_params)
	{
	  throw std::invalid_argument("Wrong number of kernel hyperparameters");
	}
      params = theta; //TODO: To make enough space. Make it more efficient.
      std::transform(theta.begin(), theta.end(), params.begin(), (double (*)(double)) exp);
    };

    vectord getHyperParameters() 
    { 
      vectord theta(params.size());
      std::transform(params.begin(), params.end(), theta.begin(), (double (*)(double)) log);
      return theta;
    };
    size_t nHyperParameters() {return n_params;};

    virtual ~AtomicKernel(){};

  protected:
    size_t n_params;
    vectord params;
  };

  //////////////////////////////////////////////////////////////////////////

  /** \brief Abstract class for isotropic kernel functors */
  class ISOkernel : public AtomicKernel
  {
  public:
    virtual ~ISOkernel(){};

  protected:
    inline double computeWeightedNorm2(const vectord &x1, const vectord &x2)
    {  
      assert(n_inputs == x1.size());
      assert(x1.size() == x2.size());
      return norm_2(x1-x2)/params(0); 
    };
  };

  /** \brief Abstract class for anisotropic kernel functors using ARD
   *   (Automatic Relevance Determination)
   */
  class ARDkernel : public AtomicKernel
  {
  public:
    virtual ~ARDkernel(){};

  protected:
    inline double computeWeightedNorm2(const vectord &x1, const vectord &x2)
    {
      assert(n_inputs == x1.size());
      assert(x1.size() == x2.size());
      assert(x1.size() == params.size());

      vectord xd = x1-x2;
      vectord r = utils::ublas_elementwise_div(xd, params);
      return norm_2(r);
    };
  };

  //@}

} //namespace bayesopt

#endif
