/**  \file criteria_hedge.hpp \brief Portfolio selection of criteria
     based on Hedge algorithm */
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

#ifndef  _CRITERIA_HEDGE_HPP_
#define  _CRITERIA_HEDGE_HPP_

#include "criteria_combined.hpp"


namespace bayesopt
{
  /**\addtogroup CriteriaFunctions*/
  //@{

  /**
   * \brief GP_Hedge model as describen in Hoffman et al. \cite Hoffman2011
   *
   * The coefficients of the bandit algorithm has been carefully selected
   * according to Shapire et al. Also, the implementation has been made to
   * avoid over or underflow.
   */
  class GP_Hedge: public CombinedCriteria
  {
  public:
    GP_Hedge();
    virtual ~GP_Hedge() {};
    void init(NonParametricProcess *proc);

    double operator() (const vectord &x) { return (*mCurrentCriterium)(x); };

    bool requireComparison(){ return true; };
    void initialCriteria();
    bool rotateCriteria();  //Returns false if NOT rotated -> end of list
    void pushResult(const vectord& prevResult);
    std::string getBestCriteria(vectord& best);

    std::string name() {return "cHedge";};
  protected:
    int update_hedge();

    vectord loss_, gain_, prob_, cumprob_;
    Criteria* mCurrentCriterium;
    std::vector<vectord> mBestLists;
    size_t mIndex;

    virtual double computeLoss(const vectord& query)
    { return mProc->prediction(query)->getMean(); }
  };

  /**
   * \brief Modification of the GP_Hedge algorithm where the bandit gains are
   * random outcomes (like Thompson sampling).
   */
  class GP_Hedge_Random: public GP_Hedge
  {
  public:
    virtual ~GP_Hedge_Random() {};
    std::string name() {return "cHedgeRandom";};
 
  protected:
    virtual double computeLoss(const vectord& query)
    { return mProc->prediction(query)->sample_query(); }
  };


  //@}

} //namespace bayesopt


#endif
