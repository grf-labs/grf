/**  \file criteria_distance.hpp \brief Cost for selecting distant points */
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

#ifndef  _CRITERIA_DISTANCE_HPP_
#define  _CRITERIA_DISTANCE_HPP_

#include "criteria_functors.hpp"

namespace bayesopt
{

  /**\addtogroup CriteriaFunctions */
  //@{

  /**
   * \brief Distance in input space. 
   * Can be combined with other critera to trade off large changes in
   * input space.
   */
  class InputDistance: public Criteria
  {
  public:
    void init(NonParametricProcess* proc)
    { 
      mProc = proc;
      mW = 1;
    };
    virtual ~InputDistance(){};
    void setParameters(const vectord &params)
    { mW = params(0); };
    size_t nParameters() {return 1;};
 
    double operator() (const vectord &x) 
    { 
      const vectord x2 = mProc->getData()->getLastSampleX();
      return mW*norm_2(x-x2);
    };
    std::string name() {return "cDistance";};
  private:
    double mW;
  };


  //@}

} //namespace bayesopt


#endif
