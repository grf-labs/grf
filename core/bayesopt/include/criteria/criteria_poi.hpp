/**  \file criteria_poi.hpp \brief Probability of improvement */
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

#ifndef  _CRITERIA_POI_HPP_
#define  _CRITERIA_POI_HPP_

#include "criteria_functors.hpp"
#include "prob_distribution.hpp"

namespace bayesopt
{

  /**\addtogroup CriteriaFunctions */
  //@{
  /// Probability of improvement criterion based on (Kushner).
  class ProbabilityOfImprovement: public Criteria
  {
  public:
    virtual ~ProbabilityOfImprovement(){};
    void init(NonParametricProcess* proc)
    { 
      mProc = proc;
      mEpsilon = 0.01;
    };
    void setParameters(const vectord &params)
    { mEpsilon = params(0); };

    size_t nParameters() {return 1;};

    inline void setEpsilon(double eps) { mEpsilon = eps; };
    double operator() (const vectord &x) 
    { 
      const double min = mProc->getValueAtMinimum();
      return mProc->prediction(x)->negativeProbabilityOfImprovement(min,
								    mEpsilon); 
    };
    std::string name() {return "cPOI";};
  private:
    double mEpsilon;
  };

  //@}

} //namespace bayesopt


#endif
