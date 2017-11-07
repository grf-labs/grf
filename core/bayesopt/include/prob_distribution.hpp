
/** \file prob_distribution.hpp 
    \brief Interface for probability models */
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


#ifndef __PROB_DISTRIBUTION_HPP__
#define __PROB_DISTRIBUTION_HPP__

#include "randgen.hpp"

namespace bayesopt
{

  class ProbabilityDistribution
  {
  public:
    ProbabilityDistribution(randEngine& eng): mtRandom(eng) {};
    virtual ~ProbabilityDistribution(){};

    /** 
     * \brief Probability density function
     * @param x query point
     * @return probability
     */
    virtual double pdf(double x) = 0;


    /** 
     * \brief Expected Improvement algorithm for minimization
     * @param min minimum value found so far
     * @param g exponent (used for annealing)
     * @return negative value of the expected improvement
     */
    virtual double negativeExpectedImprovement(double min, size_t g) = 0;

    /** 
     * \brief Lower confindence bound. Can be seen as the inverse of the Upper 
     * confidence bound
     * @param beta std coefficient (used for annealing)
     * @return value of the lower confidence bound
     */
    virtual double lowerConfidenceBound(double beta = 1) = 0;

    /** 
     * \brief Probability of improvement algorithm for minimization
     * @param min minimum value found so far
     * @param epsilon minimum improvement margin
     * 
     * @return negative value of the probability of improvement
     */
    virtual double negativeProbabilityOfImprovement(double yMin,
						    double epsilon) = 0;

    /** 
     * \brief Sample outcome acording to the marginal distribution at 
     * the query point.
     *
     * @param eng boost.random engine
     * @return outcome
     */
    virtual double sample_query() = 0;

    virtual double getMean() = 0;
    virtual double getStd() = 0;

  protected:
    randEngine& mtRandom;
  };

} //namespace bayesopt

#endif

