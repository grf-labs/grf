/** \file kernel_combined.hpp 
    \brief Kernel functions that combine other kernels */
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

#ifndef  _KERNEL_COMBINED_HPP_
#define  _KERNEL_COMBINED_HPP_

#include <boost/numeric/ublas/vector_proxy.hpp>
#include "kernel_functors.hpp"

namespace bayesopt
{
  
  /**\addtogroup KernelFunctions */
  //@{

  /** \brief Abstract class for combined kernel.
   *  It allows combinations of other kernels (addition, product, etc.)
   */
  class CombinedKernel : public Kernel
  {
  public:
    virtual void init(size_t input_dim, Kernel* left, Kernel* right)
    {
      n_inputs = input_dim;
      this->left = left;
      this->right = right;
    };

    void setHyperParameters(const vectord &theta) 
    {
      using boost::numeric::ublas::subrange;

      size_t n_lhs = left->nHyperParameters();
      size_t n_rhs = right->nHyperParameters();
      if (theta.size() != n_lhs + n_rhs)
	{
	  FILE_LOG(logERROR) << "Wrong number of kernel hyperparameters"; 
	  throw std::invalid_argument("Wrong number of kernel hyperparameters");
	}
      left->setHyperParameters(subrange(theta,0,n_lhs));
      right->setHyperParameters(subrange(theta,n_lhs,n_lhs+n_rhs));
    };

    vectord getHyperParameters() 
    {
      using boost::numeric::ublas::subrange;

      size_t n_lhs = left->nHyperParameters();
      size_t n_rhs = right->nHyperParameters();
      vectord par(n_lhs + n_rhs);
      subrange(par,0,n_lhs) = left->getHyperParameters();
      subrange(par,n_lhs,n_lhs+n_rhs) = right->getHyperParameters();
      return par;
    };

    size_t nHyperParameters() 
    {
      size_t n_lhs = left->nHyperParameters();
      size_t n_rhs = right->nHyperParameters();
      return n_lhs + n_rhs;
    };

    virtual ~CombinedKernel()
    {
      delete left;
      delete right;
    };

  protected:
    Kernel* left;
    Kernel* right;
  };

  //@}

} //namespace bayesopt

#endif
