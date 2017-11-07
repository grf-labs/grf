/** \file mean_combined.hpp 
    \brief Parametric functions that combine other functions */
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

#ifndef  _MEAN_COMBINED_HPP_
#define  _MEAN_COMBINED_HPP_

#include <boost/numeric/ublas/vector_proxy.hpp>
#include "mean_functors.hpp"

namespace bayesopt
{

  /**\addtogroup ParametricFunctions
   * @{
   */

  /** \brief Abstract class for combined functions.
   *  It allows combinations of other functions (addition, product, etc.)
   */
  class CombinedFunction : public ParametricFunction
  {
  public:
    virtual int init(size_t input_dim, 
		     ParametricFunction* left, 
		     ParametricFunction* right)
    {
      n_inputs = input_dim;
      this->left = left;
      this->right = right;
      return 0;
    };
    void setParameters(const vectord &theta) 
    {
      using boost::numeric::ublas::subrange;

      size_t n_lhs = left->nParameters();
      size_t n_rhs = right->nParameters();
      if (theta.size() != n_lhs + n_rhs)
	{
	  throw std::invalid_argument("Wrong number of mean function parameters"); 
	}

      left->setParameters(subrange(theta,0,n_lhs));
      right->setParameters(subrange(theta,n_lhs,n_lhs+n_rhs));
    };

    vectord getParameters() 
    {
      using boost::numeric::ublas::subrange;

      size_t n_lhs = left->nParameters();
      size_t n_rhs = right->nParameters();
      vectord par(n_lhs + n_rhs);
      subrange(par,0,n_lhs) = left->getParameters();
      subrange(par,n_lhs,n_lhs+n_rhs) = right->getParameters();
      return par;
    };

    size_t nParameters() 
    {
      size_t n_lhs = left->nParameters();
      size_t n_rhs = right->nParameters();
      return n_lhs + n_rhs;
    };

    virtual ~CombinedFunction()
    {
      delete left;
      delete right;
    };

  protected:
    ParametricFunction* left;
    ParametricFunction* right;
  };


  /** \brief Sum of two kernels */
  class SumFunction: public CombinedFunction
  {
  public:
    double getMean(const vectord &x)
    {
      return left->getMean(x) + right->getMean(x);
    };

    vectord getFeatures(const vectord &x)
    {
      using boost::numeric::ublas::subrange;

      size_t n_lhf = left->nFeatures();
      size_t n_rhf = right->nFeatures();
      vectord feat(n_lhf + n_rhf);
      subrange(feat,0,n_lhf) = left->getFeatures(x);
      subrange(feat,n_lhf,n_lhf+n_rhf) = right->getFeatures(x);
      return feat;
    };

    size_t nFeatures() 
    {
      size_t n_lhf = left->nFeatures();
      size_t n_rhf = right->nFeatures();
      return n_lhf + n_rhf;
    };

  };

  //@}

} //namespace bayesopt



#endif
