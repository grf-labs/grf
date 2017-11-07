/*
-----------------------------------------------------------------------------
   This file is part of BayesOptimization, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2015 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   BayesOptimization is free software: you can redistribute it and/or modify
   it under the terms of the GNU Affero General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   BayesOptimization is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Affero General Public License for more details.

   You should have received a copy of the GNU Affero General Public License
   along with BayesOptimization.  If not, see <http://www.gnu.org/licenses/>.
-----------------------------------------------------------------------------
*/

#include "bayesopt/bayesopt.h"
#include "bayesopt/bayesopt.hpp"      

#include "log.hpp"
#include "ublas_extra.hpp"
#include "specialtypes.hpp"

static const int BAYESOPT_FAILURE = -1; /* generic failure code */
static const int BAYESOPT_INVALID_ARGS = -2;
static const int BAYESOPT_OUT_OF_MEMORY = -3;
static const int BAYESOPT_RUNTIME_ERROR = -4;

/**
 * \brief Version of ContinuousModel for the C wrapper
 */
class CContinuousModel: public bayesopt::ContinuousModel 
{
 public:

  CContinuousModel(size_t dim, bopt_params params):
    ContinuousModel(dim,params)  {}; 

  virtual ~CContinuousModel(){};

  double evaluateSample( const vectord &Xi ) 
  {
    int n = static_cast<int>(Xi.size());
    return  mF(n,&Xi[0],NULL,mOtherData);
  };

  void set_eval_funct(eval_func f)
  {  mF = f; }


  void save_other_data(void* other_data)
  {  mOtherData = other_data; }
 
protected:
  void* mOtherData;
  eval_func mF;
};

/**
 * \brief Version of DiscreteModel for the C wrapper
 */
class CDiscreteModel: public bayesopt::DiscreteModel
{
 public:

  CDiscreteModel(const vecOfvec &validX, bopt_params params):
    DiscreteModel(validX, params)
  {}; 

  CDiscreteModel(const vectori &categories, bopt_params params):
    DiscreteModel(categories, params)
  {}; 

  double evaluateSample( const vectord &Xi ) 
  {
    int n = static_cast<int>(Xi.size());
    return  mF(n,&Xi[0],NULL,mOtherData);
  };

  void set_eval_funct(eval_func f)
  {  mF = f; }


  void save_other_data(void* other_data)
  {  mOtherData = other_data; }
 
protected:
  void* mOtherData;
  eval_func mF;
};

int bayes_optimization(int nDim, eval_func f, void* f_data,
		       const double *lb, const double *ub,
		       double *x, double *minf, bopt_params parameters)
{
  vectord result(nDim);

  vectord lowerBound = bayesopt::utils::array2vector(lb,nDim); 
  vectord upperBound = bayesopt::utils::array2vector(ub,nDim); 

  try 
    {
      CContinuousModel optimizer(nDim, parameters);

      optimizer.set_eval_funct(f);
      optimizer.save_other_data(f_data);
      optimizer.setBoundingBox(lowerBound,upperBound);

      optimizer.optimize(result);
      std::copy(result.begin(), result.end(), x);

      *minf = optimizer.getValueAtMinimum();
    }
  catch (std::bad_alloc& e)
    {
      FILE_LOG(logERROR) << e.what(); 
      return  BAYESOPT_OUT_OF_MEMORY; 
    }
  catch (std::invalid_argument& e)
    { 
      FILE_LOG(logERROR) << e.what(); 
      return BAYESOPT_INVALID_ARGS; 
    }
  catch (std::runtime_error& e)
    { 
      FILE_LOG(logERROR) << e.what(); 
      return BAYESOPT_RUNTIME_ERROR;
    }
  catch (...)
    { 
      FILE_LOG(logERROR) << "Unknown error";
      return BAYESOPT_FAILURE; 
    }
  return 0; /* everything ok*/
};

int bayes_optimization_disc(int nDim, eval_func f, void* f_data,
			    double *valid_x, size_t n_points,
			    double *x, double *minf, bopt_params parameters)
{
  vectord result(nDim);
  vectord input(nDim);
  vecOfvec xSet;

  for(size_t i = 0; i<n_points;++i)
    {
      for(int j = 0; j<nDim; ++j)
	{
	 input(j) = valid_x[i*nDim+j]; 
	}
      xSet.push_back(input);
    }

  if(parameters.n_init_samples > n_points)
    {
      parameters.n_init_samples = n_points;
      parameters.n_iterations = 0;
    }

  try
    {
      CDiscreteModel optimizer(xSet,parameters);
      
      optimizer.set_eval_funct(f);
      optimizer.save_other_data(f_data);
      optimizer.optimize(result);

      std::copy(result.begin(), result.end(), x);

      *minf = optimizer.getValueAtMinimum();
    }
  catch (std::bad_alloc& e)
    {
      FILE_LOG(logERROR) << e.what(); 
      return  BAYESOPT_OUT_OF_MEMORY; 
    }
  catch (std::invalid_argument& e)
    { 
      FILE_LOG(logERROR) << e.what(); 
      return BAYESOPT_INVALID_ARGS; 
    }
  catch (std::runtime_error& e)
    { 
      FILE_LOG(logERROR) << e.what(); 
      return BAYESOPT_RUNTIME_ERROR;
    }
  catch (...)
    { 
      FILE_LOG(logERROR) << "Unknown error";
      return BAYESOPT_FAILURE; 
    }

  return 0; /* everything ok*/
}


int bayes_optimization_categorical(int nDim, eval_func f, void* f_data,
				   int *categories, double *x, 
				   double *minf, bopt_params parameters)
{
  vectord result(nDim);
  vectori cat(nDim);

  std::copy(categories,categories+nDim,cat.begin());

  try
    {
      CDiscreteModel optimizer(cat,parameters);
      
      optimizer.set_eval_funct(f);
      optimizer.save_other_data(f_data);
      optimizer.optimize(result);

      std::copy(result.begin(), result.end(), x);

      *minf = optimizer.getValueAtMinimum();
    }
  catch (std::bad_alloc& e)
    {
      FILE_LOG(logERROR) << e.what(); 
      return  BAYESOPT_OUT_OF_MEMORY; 
    }
  catch (std::invalid_argument& e)
    { 
      FILE_LOG(logERROR) << e.what(); 
      return BAYESOPT_INVALID_ARGS; 
    }
  catch (std::runtime_error& e)
    { 
      FILE_LOG(logERROR) << e.what(); 
      return BAYESOPT_RUNTIME_ERROR;
    }
  catch (...)
    { 
      FILE_LOG(logERROR) << "Unknown error";
      return BAYESOPT_FAILURE; 
    }

  return 0; /* everything ok*/
}
