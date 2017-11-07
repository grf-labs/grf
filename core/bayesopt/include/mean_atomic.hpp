/** \file mean_atomic.hpp \brief Atomic (simple) parametric functions */
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

#ifndef  _MEAN_ATOMIC_HPP_
#define  _MEAN_ATOMIC_HPP__

#include <boost/numeric/ublas/vector_proxy.hpp>
#include "mean_functors.hpp"

namespace bayesopt
{

  /**\addtogroup ParametricFunctions
   * @{
   */

  /** \brief Abstract class for an atomic kernel */
  class AtomicFunction : public ParametricFunction
  {
  public:
    virtual int init(size_t input_dim)
    {
      n_inputs = input_dim;
      return 0;
    };
    void setParameters(const vectord &theta) 
    {
      if(theta.size() != n_params)
	{
	  throw std::invalid_argument("Wrong number of mean function parameters"); 
	}
   
      mParameters = theta;
    };
    vectord getParameters() {return mParameters;};
    size_t nParameters() {return n_params;};
    size_t nFeatures() {return n_features;};

    virtual ~AtomicFunction(){};

  protected:
    size_t n_params;
    size_t n_features;
    vectord mParameters;
  };


  /** \brief Constant zero function */
  class ZeroFunction: public AtomicFunction
  {
  public:
    int init(size_t input_dim)
    {
      n_inputs = input_dim;
      n_params = 1;
      n_features = 1;
      return 0;
    };
    double getMean (const vectord& x) { return 0.0; };
    vectord getFeatures(const vectord& x) { return zvectord(1); };  
  };

  /** \brief Constant one function */
  class OneFunction: public AtomicFunction
  {
  public:
    int init(size_t input_dim)
    {
      n_inputs = input_dim;
      n_params = 1;
      n_features = 1;
      return 0;
    };
    double getMean (const vectord& x) { return 1.0; };
    vectord getFeatures(const vectord& x) { return svectord(1,1.0); };  
  };


  /** \brief Constant function. 
      The first parameter indicates the constant value. */
  class ConstantFunction: public AtomicFunction
  {
  public:
    int init(size_t input_dim)
    {
      n_inputs = input_dim;
      n_params = 1;
      n_features = 1;
      return 0;
    };
    double getMean (const vectord& x) { return mParameters(0); };
    vectord getFeatures(const vectord& x) { return svectord(1,1.0); };  
  };


  /** \brief Linear combination function. 
      Each parameter indicates the coefficient of each dimension. */
  class LinearFunction: public AtomicFunction
  {
  public:
    int init(size_t input_dim)
    {
      n_inputs = input_dim;
      n_params = input_dim;
      n_features = input_dim;
      return 0;
    };
    double getMean (const vectord& x)
    { return boost::numeric::ublas::inner_prod(x,mParameters);  };
    vectord getFeatures(const vectord& x) { return x; };  
  };


  /** \brief Linear combination plus constant function. 
      The first parameter indicates the constant value. */
  class LinearPlusConstantFunction: public AtomicFunction
  {
  public:
    int init(size_t input_dim)
    {
      n_inputs = input_dim;
      n_params = input_dim + 1;
      n_features = input_dim + 1;
      return 0;
    };
    void setParameters(const vectord& params)
    { 
      if(params.size() != n_params)
	{
	  throw std::invalid_argument("Wrong number of mean function parameters"); 
	}

      mConstParam = params(0);
      mParameters = boost::numeric::ublas::project(params, 
						   boost::numeric::ublas::range(1, params.size())); 
    };
  
    double getMean (const vectord& x)
    { return boost::numeric::ublas::inner_prod(x,mParameters) + mConstParam;  };

    vectord getFeatures(const vectord& x) 
    {
      using boost::numeric::ublas::range;
      using boost::numeric::ublas::project;
      vectord res(x.size()+1);
      res(0) = 1;
      project(res,range(1,res.size())) = x;
      return res; 
    };  

  protected:
    double mConstParam;
  };

  //@}

} //namespace bayesopt

#endif
