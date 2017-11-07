/**  \file criteria_a_opt.hpp \brief A-optimality (uncertainty) based criteria */
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

#ifndef  _CRITERIA_A_OPT_HPP_
#define  _CRITERIA_A_OPT_HPP_

#include "criteria_functors.hpp"

namespace bayesopt
{

  /**\addtogroup CriteriaFunctions */
  //@{

  /**
   * \brief Greedy A-Optimality criterion.  
   * Used for learning the function, not to minimize. Some authors
   * name it I-optimality because it minimizes the error on the
   * prediction, not on the parameters.
   */
  class GreedyAOptimality: public Criteria
  {
  public:
    virtual ~GreedyAOptimality(){};
    void setParameters(const vectord &params) {};
    size_t nParameters() {return 0;};
    double operator() (const vectord &x) 
    { return -mProc->prediction(x)->getStd(); };
    std::string name() {return "cAopt";};
  };

  //@}

} //namespace bayesopt


#endif
