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

#ifndef  _CRITERIA_LCB_HPP_
#define  _CRITERIA_LCB_HPP_

#include "criteria_functors.hpp"

namespace bayesopt
{

  /**\addtogroup CriteriaFunctions  */
  //@{

  /// Lower (upper) confidence bound criterion by [Cox and John, 1992].
  class LowerConfidenceBound: public Criteria
  {
  public:
    virtual ~LowerConfidenceBound(){};
    void init(NonParametricProcess* proc)
    { 
      mProc = proc;
      mBeta = 1.0;
    };
    void setParameters(const vectord &params)
    { mBeta = params(0); };

    size_t nParameters() {return 1;};

    double operator() (const vectord &x) 
    { 
      return mProc->prediction(x)->lowerConfidenceBound(mBeta); 
    };
    std::string name() {return "cLCB";};
  private:
    double mBeta;
  };


  /// Lower (upper) confidence bound using Srinivas annealing \cite Srinivas10
  class AnnealedLowerConfindenceBound: public Criteria
  {
  public:
    virtual ~AnnealedLowerConfindenceBound(){};
    void init(NonParametricProcess* proc)
    { 
      mProc = proc;
      reset();
    };

    void setParameters(const vectord &params)
    { mCoef = params(0); };

    size_t nParameters() {return 1;};
    void reset() { nCalls = 1; mCoef = 5.0;};
    double operator() (const vectord &x) 
    {
      size_t nDims = x.size();
    
      double beta = sqrt(2*log(static_cast<double>(nCalls*nCalls))*(nDims+1) 
			 + log(static_cast<double>(nDims))*nDims*mCoef);
      ProbabilityDistribution* d_ = mProc->prediction(x);
      return d_->lowerConfidenceBound(beta); 
    };
    void update(const vectord &x) { ++nCalls; }

    std::string name() {return "cLCBa";};
  private:
    double mCoef;
    unsigned int nCalls;
  };

  //@}

} //namespace bayesopt


#endif
