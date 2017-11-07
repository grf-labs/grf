/**  \file posterior_fixed.hpp \brief Posterior model based on fixed kernel parameters */
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


#ifndef  _POSTERIOR_FIXED_HPP_
#define  _POSTERIOR_FIXED_HPP_

#include <boost/scoped_ptr.hpp>
#include "inneroptimization.hpp"
#include "criteria_functors.hpp"
#include "posteriormodel.hpp"


namespace bayesopt {

  /** \addtogroup BayesOpt
   *  \brief Main module for Bayesian optimization
   */
  /*@{*/

 /**
   * \brief Bayesian optimization using different non-parametric 
   * processes as distributions over surrogate functions. 
   */
  class PosteriorFixed: public PosteriorModel
  {
  public:
    /** 
     * Constructor
     * @param params set of parameters (see parameters.hpp)
     */
    PosteriorFixed(size_t dim, Parameters params, randEngine& eng);

    /** 
     * Default destructor
     */
    virtual ~PosteriorFixed();

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
    PosteriorFixed();

    void setSurrogateModel(randEngine& eng);    
    void setCriteria(randEngine& eng);

  private:  // Members
    boost::scoped_ptr<NonParametricProcess> mGP; ///< Pointer to surrogate model
    boost::scoped_ptr<Criteria> mCrit;                   ///< Metacriteria model
  };

  /**@}*/

  inline void PosteriorFixed::updateHyperParameters()
  { /* Fixed model. Does not require updates */ };

  inline void PosteriorFixed::fitSurrogateModel()
  { mGP->fitSurrogateModel(); };

  inline void PosteriorFixed::updateSurrogateModel()
  { mGP->updateSurrogateModel(); };

  inline double PosteriorFixed::evaluateCriteria(const vectord& query)
  { return (*mCrit)(query); };

  inline void PosteriorFixed::updateCriteria(const vectord& query)
  { return mCrit->update(query); };

  inline bool PosteriorFixed::criteriaRequiresComparison()
  {return mCrit->requireComparison(); };
    
  inline void PosteriorFixed::setFirstCriterium()
  { mCrit->initialCriteria(); };

  inline bool PosteriorFixed::setNextCriterium(const vectord& prevResult)
  { 
    mCrit->pushResult(prevResult);
    return mCrit->rotateCriteria(); 
  };

  inline std::string PosteriorFixed::getBestCriteria(vectord& best)
  { return mCrit->getBestCriteria(best); };

  inline  ProbabilityDistribution* PosteriorFixed::getPrediction(const vectord& query)
  { return mGP->prediction(query); };


} //namespace bayesopt


#endif
