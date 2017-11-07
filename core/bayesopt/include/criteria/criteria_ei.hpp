/**  \file criteria_ei.hpp \brief Expected improvement based criteria */
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

#ifndef  _CRITERIA_EI_HPP_
#define  _CRITERIA_EI_HPP_

#include "criteria_functors.hpp"
#include "prob_distribution.hpp"

namespace bayesopt
{

  /**\addtogroup CriteriaFunctions  */
  //@{

  /// Expected improvement criterion by Mockus \cite Mockus78
  class ExpectedImprovement: public Criteria
  {
  public:
    virtual ~ExpectedImprovement(){};
    void init(NonParametricProcess* proc)
    { 
      mProc = proc;
      mExp = 1;
    };

    void setParameters(const vectord &params)
    { mExp = static_cast<size_t>(params(0)); };

    size_t nParameters() {return 1;};

    double operator() (const vectord &x) 
    { 
      const double min = mProc->getValueAtMinimum();
      return mProc->prediction(x)->negativeExpectedImprovement(min,mExp); 
    };

    std::string name() {return "cEI";};

  private:
    size_t mExp;
  };

  /// Expected improvement criterion modification by Lizotte
  class BiasedExpectedImprovement: public Criteria
  {
  public:
    virtual ~BiasedExpectedImprovement(){};
    void init(NonParametricProcess* proc)
    { 
      mProc = proc;
      mBias = 0.01;
      mExp = 1;
    };

    void setParameters(const vectord &params)
    {
      mExp = static_cast<size_t>(params(0));
      mBias = params(1);
    };

    size_t nParameters() {return 2;};

    double operator() (const vectord &x) 
    { 
      const double sigma = mProc->getSignalVariance();
      const double min = mProc->getValueAtMinimum() - mBias/sigma;
      return mProc->prediction(x)->negativeExpectedImprovement(min,mExp); 
    };
    std::string name() {return "cBEI";};
  private:
    double mBias;
    size_t mExp;
  };


  /// Expected improvement criterion using Schonlau annealing. \cite Schonlau98
  class AnnealedExpectedImprovement: public Criteria
  {
  public:
    virtual ~AnnealedExpectedImprovement(){};
    void init(NonParametricProcess* proc)
    { 
      mProc = proc;
      reset();
    };

    void setParameters(const vectord &params)
    { mExp = static_cast<size_t>(params(0)); };

    size_t nParameters() {return 1;};
    void reset() { nCalls = 1; mExp = 10;};
    double operator() (const vectord &x) 
    {
      ProbabilityDistribution* d_ = mProc->prediction(x);
      const double min = mProc->getValueAtMinimum();
      return d_->negativeExpectedImprovement(min,mExp); 
    };
    void update(const vectord &x) 
    {
      ++nCalls;
      if (nCalls % 10)
	mExp = static_cast<size_t>(ceil(mExp/2.0));
    }
    std::string name() {return "cEIa";};
  private:
    size_t mExp;
    unsigned int nCalls;
  };



  //@}

} //namespace bayesopt


#endif
