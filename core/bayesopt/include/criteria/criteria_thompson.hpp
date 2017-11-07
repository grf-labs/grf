/**  \file criteria_thompson.hpp \brief Thompson and optimistic
     sampling criteria */
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

#ifndef  _CRITERIA_THOMPSON_HPP_
#define  _CRITERIA_THOMPSON_HPP_

#include "criteria_functors.hpp"

namespace bayesopt
{

  /**\addtogroup CriteriaFunctions */
  //@{

  /// Thompson sampling. Picks a random sample of the surrogate model.
  class ThompsonSampling: public Criteria
  {
  public:
    ThompsonSampling() {};
    virtual ~ThompsonSampling(){};
    void setParameters(const vectord &params) { };
    size_t nParameters() {return 0;};
    double operator() (const vectord &x) 
    {
      ProbabilityDistribution* d_ = mProc->prediction(x);
      return d_->sample_query();
    };
    std::string name() {return "cThompsonSampling";};
  };

  /**
   * \brief Optimistic sampling. 
   * A simple variation of Thompson sampling that picks only samples
   * that are better than the best outcome so far.
   */
  class OptimisticSampling: public Criteria
  {
  public:
    OptimisticSampling() {};
    virtual ~OptimisticSampling(){};
    void setParameters(const vectord &params) {};
    size_t nParameters() {return 0;};
    double operator() (const vectord &x)  
    {
      ProbabilityDistribution* d_ = mProc->prediction(x);
      const double yStar = d_->sample_query();
      const double yPred = d_->getMean();
      return (std::min)(yPred,yStar);
    };
    std::string name() {return "cOptimisticSampling";};
  };


  //@}

} //namespace bayesopt


#endif
