/**  \file criteria_functors.hpp \brief Abstract and factory modules
     for criteria */
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

#ifndef  _CRITERIA_FUNCTORS_HPP_
#define  _CRITERIA_FUNCTORS_HPP_

#include <map>
#include "nonparametricprocess.hpp"

namespace bayesopt
{

  /**\addtogroup CriteriaFunctions
   * \brief Set of criterium functions to select the next point during
   * optimization/exploration
   */
  //@{

  /**
   * \brief Abstract interface for criteria functors.
   */
  class Criteria
  {
  public:
    virtual ~Criteria() {};
    virtual void init(NonParametricProcess *proc) { mProc = proc; };

    double evaluate(const vectord &x)  {return (*this)(x);}
    virtual double operator() (const vectord &x)  = 0;

    virtual std::string name() = 0;
    virtual void setParameters(const vectord &params) = 0;
    virtual size_t nParameters() = 0;

    //Dummy functions. Not all criteria support these methods.
    virtual void reset() { assert(false); };
    void setRandomEngine(randEngine& eng){ mtRandom = &eng; }

    // In general, most criteria does not support comparisons!
    // TODO: Consider throwing exception when incorrect calls
    virtual void pushCriteria(Criteria* crit){};
    virtual bool requireComparison(){ return false; };
    virtual void initialCriteria(){};
    virtual void update(const vectord &x){};
    virtual bool rotateCriteria(){return false;};
    virtual void pushResult(const vectord& prevResult){};
    virtual std::string getBestCriteria(vectord& best)
    { assert(false); return name(); };

  protected:
    NonParametricProcess *mProc;
    randEngine* mtRandom;
  };


  /** 
   * \brief Factory model for criterion functions
   * This factory is based on the libgp library by Manuel Blum
   *      https://bitbucket.org/mblum/libgp
   * which follows the squeme of GPML by Rasmussen and Nickisch
   *     http://www.gaussianprocess.org/gpml/code/matlab/doc/
   */
  class CriteriaFactory
  {
  public:
    CriteriaFactory ();
    virtual ~CriteriaFactory () {};
  
    //Criteria* create(criterium_name name, NonParametricProcess* proc);
    Criteria* create(std::string name, NonParametricProcess* proc);
    
  private:
    typedef Criteria* (*create_func_definition)();
    std::map<std::string , CriteriaFactory::create_func_definition> registry;
  };


  //@}

} //namespace bayesopt


#endif
