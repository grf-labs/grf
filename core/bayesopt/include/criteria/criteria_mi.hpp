/**  \file criteria_lcb.hpp \brief Lower confidence bound based criteria */
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

#ifndef  _CRITERIA_MI_HPP_
#define  _CRITERIA_MI_HPP_

#include "criteria_functors.hpp"

namespace bayesopt
{

  /**\addtogroup CriteriaFunctions  */
  //@{

  /// Mutual Information bound criterion b [Contal et al., 2014].
  class MutualInformation: public Criteria
  {
  public:
    virtual ~MutualInformation(){};
    void init(NonParametricProcess* proc)
    { 
      mProc = proc;
      mSqAlpha = sqrt(std::log(2/1e-6));  // See [Contal et al., 2014].
      mGamma = 0.0;
    };
    void setParameters(const vectord &params)
    { mSqAlpha = sqrt(params(0)); };

    size_t nParameters() {return 1;};

    double operator() (const vectord &x) 
    { 
      ProbabilityDistribution* d = mProc->prediction(x);
      double mu = d->getMean();
      double sigma2 = d->getStd() * d->getStd();
      return mu + mSqAlpha * (sqrt(sigma2+mGamma) - sqrt(mGamma));
    };
    void update(const vectord &x)
    {
      ProbabilityDistribution* d = mProc->prediction(x);
      double mu = d->getMean();
      double sigma2 = d->getStd() * d->getStd();
      mGamma += sigma2; 
    }
    std::string name() {return "cMI";};
  private:
    double mSqAlpha;
    double mGamma;
  };

  //@}

} //namespace bayesopt


#endif
