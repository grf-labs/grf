/** \file gauss_distribution.hpp 
    \brief Gaussian probability distribution */
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


#ifndef __GAUSS_DISTRIBUTION_HPP__
#define __GAUSS_DISTRIBUTION_HPP__

#include <boost/math/distributions/normal.hpp> 
#include "prob_distribution.hpp"

namespace bayesopt
{

  class GaussianDistribution: public ProbabilityDistribution
  {
  public:
    GaussianDistribution(randEngine& eng);
    virtual ~GaussianDistribution();

    /** 
     * \brief Sets the mean and std of the distribution
     */
    void setMeanAndStd(double mean, double std)
    { mean_ = mean; std_ = std; };

    /** 
     * \brief Probability density function
     * @param x query point
     * @return probability
     */
    double pdf(double x) 
    {
      x = (x - mean_) / std_;
      return boost::math::pdf(d_,x); 
    };

    /** 
     * \brief Expected Improvement algorithm for minimization
     * @param min  minimum value found
     * @param g exponent (used for annealing)
     *
     * @return negative value of the expected improvement
     */
    double negativeExpectedImprovement(double min, size_t g);

    /** 
     * \brief Lower confindence bound. Can be seen as the inverse of the Upper 
     * confidence bound
     * @param beta std coefficient (used for annealing)
     * @return value of the lower confidence bound
     */
    double lowerConfidenceBound(double beta);

    /** 
     * Probability of improvement algorithm for minimization
     * @param min  minimum value found
     * @param epsilon minimum improvement margin
     * 
     * @return negative value of the probability of improvement
     */
    double negativeProbabilityOfImprovement(double min,
					    double epsilon);

    /** 
     * Sample outcome acording to the marginal distribution at the query point.
     * @return outcome
     */
    double sample_query();

    double getMean() { return mean_; };
    double getStd()  { return std_; };


  private:
    boost::math::normal d_;
    double mean_;
    double std_;
  };

} //namespace bayesopt

#endif
