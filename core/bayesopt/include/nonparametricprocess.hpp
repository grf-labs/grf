/** \file nonparametricprocess.hpp 
    \brief Abstract module for a Bayesian regressor */
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


#ifndef __BAYESIANREGRESSOR_HPP__
#define __BAYESIANREGRESSOR_HPP__

#include "dataset.hpp"
#include "prob_distribution.hpp"
#include "mean_functors.hpp"
#include "optimizable.hpp"

namespace bayesopt
{

  /** \addtogroup NonParametricProcesses */
  /**@{*/

  /**
   * \brief Abstract class to implement Bayesian regressors
   */
  class NonParametricProcess: public RBOptimizable
  {
  public:
    NonParametricProcess(size_t dim, Parameters parameters, 
			 const Dataset& data, 
			 MeanModel& mean,
			 randEngine& eng);

    virtual ~NonParametricProcess();

    /** 
     * \brief Factory model generator for surrogate models
     * @param parameters (process name, noise, priors, etc.)
     * @return pointer to the corresponding derivate class (surrogate model)
     */
    static NonParametricProcess* create(size_t dim, Parameters parameters,
					const Dataset& data, 			 
					MeanModel& mean, randEngine& eng);

    /** 
     * \brief Function that returns the prediction of the GP for a query point
     * in the hypercube [0,1].
     * 
     * @param query in the hypercube [0,1] to evaluate the Gaussian process
     * @return pointer to the probability distribution.
     */	
    virtual ProbabilityDistribution* prediction(const vectord &query) = 0;
		 		 
    /** 
     * \brief Computes the initial surrogate model and updates the
     * kernel parameters estimation. 
     *
     * This function requires to recompute all covariance matrixes,
     * inverses, etc.  Use it with precaution.
     */
    virtual void fitSurrogateModel() = 0;

    /** 
     * \brief Sequential update of the surrogate model by adding a new
     * row to the Kernel matrix, more precisely, to its Cholesky
     * decomposition. 
     * 
     * It assumes that the kernel hyperparemeters do not change.
     */   
    virtual void updateSurrogateModel() = 0;


    // Getters and setters
    double getValueAtMinimum();
    const Dataset* getData();
    double getSignalVariance();
    
    virtual size_t nHyperParameters() = 0;
    virtual vectord getHyperParameters() = 0;
    virtual void setHyperParameters(const vectord& theta) = 0;


  protected:
    const Dataset& mData;  
    double mSigma;                                   //!< Signal variance
    size_t dim_;
    MeanModel& mMean;
  };

  //////////////////////////////////////////////////////////////////////////////
  //// Inlines

  inline double NonParametricProcess::getValueAtMinimum() 
  { return mData.getValueAtMinimum(); };


  inline const Dataset* NonParametricProcess::getData()
  {return &mData;}

  inline double NonParametricProcess::getSignalVariance() 
  { return mSigma; };

  /**@}*/
  
} //namespace bayesopt

#endif
