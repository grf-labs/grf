/** \file gaussian_process.hpp 
    \brief Standard zero mean gaussian process with noisy observations */
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

#ifndef  _GAUSSIAN_PROCESS_HPP_
#define  _GAUSSIAN_PROCESS_HPP_

#include "gauss_distribution.hpp"
#include "conditionalbayesprocess.hpp"


namespace bayesopt
{
  
  /** \addtogroup NonParametricProcesses */
  /**@{*/

  /**
   * \brief Standard zero mean gaussian process with noisy observations.
   */
  class GaussianProcess: public ConditionalBayesProcess
  {
  public:
    GaussianProcess(size_t dim, Parameters params, const Dataset& data, 
		    MeanModel& mean, randEngine& eng);
    virtual ~GaussianProcess();

    /** 
     * \brief Function that returns the prediction of the GP for a query point
     * in the hypercube [0,1].
     * 
     * @param query in the hypercube [0,1] to evaluate the Gaussian process
     * @return pointer to the probability distribution.
     */	
    ProbabilityDistribution* prediction(const vectord &query);

  private:

    /** 
     * \brief Computes the negative log likelihood of the data for all
     * the parameters.
     * @return value negative log likelihood
     */
    double negativeTotalLogLikelihood();

    /** 
     * \brief Computes the negative log likelihood of the data.
     * 
     * \f[ \log p(y|x,\theta,f) \propto  y^T (K+\sigma I)^{-1} y + 
     *                             \log|K+\sigma I|
     * \f]
     *
     * @return value negative log likelihood
     */
    double negativeLogLikelihood();

    /** Precompute some values of the prediction that do not depends
     *	on the query
     */
    void precomputePrediction();

  private:
    vectord mAlphaV;              ///< Precomputed L\y
    GaussianDistribution* d_;     ///< Pointer to distribution function
  };

  /**@}*/

} //namespace bayesopt
 

#endif
