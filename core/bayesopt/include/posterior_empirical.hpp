/**  \file empiricalbayes.hpp \brief Empirical Bayesian estimator */
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


#ifndef  _EMPIRICALBAYES_HPP_
#define  _EMPIRICALBAYES_HPP_

#include <boost/scoped_ptr.hpp>
#include "criteria_functors.hpp"
#include "posteriormodel.hpp"

namespace bayesopt {

  /** \addtogroup BayesOpt
   *  \brief Main module for Bayesian optimization
   */
  /*@{*/

  class NLOPT_Optimization;

 /**
   * \brief Bayesian optimization using different non-parametric 
   * processes as distributions over surrogate functions. 
   */
  class EmpiricalBayes: public PosteriorModel
  {
  public:
    /** 
     * Constructor
     * @param params set of parameters (see parameters.hpp)
     */
    EmpiricalBayes(size_t dim, Parameters params, randEngine& eng);

    /** 
     * Default destructor
     */
    virtual ~EmpiricalBayes();

    void updateHyperParameters();
    void fitSurrogateModel();
    void updateSurrogateModel();

    double evaluateCriteria(const vectord& query);
    void updateCriteria(const vectord& query);

    bool criteriaRequiresComparison();
    void setFirstCriterium();
    bool setNextCriterium(const vectord& prevResult);
    std::string getBestCriteria(vectord& best);

    ProbabilityDistribution* getPrediction(const vectord& query);

  private:
    EmpiricalBayes();

    void setSurrogateModel(randEngine& eng);    
    void setCriteria(randEngine& eng);

  private:  // Members
    boost::scoped_ptr<NonParametricProcess> mGP; ///< Pointer to surrogate model
    boost::scoped_ptr<Criteria> mCrit;                   ///< Metacriteria model

    boost::scoped_ptr<NLOPT_Optimization> kOptimizer;
  };

  /**@}*/

  inline void EmpiricalBayes::fitSurrogateModel()
  { mGP->fitSurrogateModel(); };

  inline void EmpiricalBayes::updateSurrogateModel()
  { mGP->updateSurrogateModel(); };

  inline double EmpiricalBayes::evaluateCriteria(const vectord& query)
  { return (*mCrit)(query); };

  inline void EmpiricalBayes::updateCriteria(const vectord& query)
  { return mCrit->update(query); };

  inline bool EmpiricalBayes::criteriaRequiresComparison()
  {return mCrit->requireComparison(); };
    
  inline void EmpiricalBayes::setFirstCriterium()
  { mCrit->initialCriteria(); };

  inline bool EmpiricalBayes::setNextCriterium(const vectord& prevResult)
  { 
    mCrit->pushResult(prevResult);
    return mCrit->rotateCriteria(); 
  };

  inline std::string EmpiricalBayes::getBestCriteria(vectord& best)
  { return mCrit->getBestCriteria(best); };

  inline  ProbabilityDistribution* EmpiricalBayes::getPrediction(const vectord& query)
  { return mGP->prediction(query); };


} //namespace bayesopt


#endif
