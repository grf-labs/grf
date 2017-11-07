/**  \file criteria_expected.hpp \brief Criterion based on the expected
     value of the function. */
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

#ifndef  _CRITERIA_EXPECTED_HPP_
#define  _CRITERIA_EXPECTED_HPP_

#include "criteria_functors.hpp"

namespace bayesopt
{

  /**\addtogroup CriteriaFunctions  */
  //@{

  /// Expected return criterion.
  class ExpectedReturn: public Criteria
  {
  public:
    virtual ~ExpectedReturn(){};
    void setParameters(const vectord &params) { };
    size_t nParameters() {return 0;};
    double operator() (const vectord &x) 
    { return mProc->prediction(x)->getMean(); };
    std::string name() {return "cExpReturn";};
  };


  //@}

} //namespace bayesopt


#endif
