/**  \file posteriormodel.hpp \brief Abstract interface for posterior model/criteria */
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


#ifndef  _POSTERIORMODEL_HPP_
#define  _POSTERIORMODEL_HPP_

#include "prob_distribution.hpp"
#include "dataset.hpp"
#include "mean_functors.hpp"

namespace bayesopt {

  /** \addtogroup BayesOpt
   *  \brief Main module for Bayesian optimization
   */
  /*@{*/

 /**
   * \brief Bayesian optimization using different non-parametric 
   * processes as distributions over surrogate functions. 
   */
  class PosteriorModel
  {
  public:
    static PosteriorModel* create(size_t dim, Parameters params, 
				  randEngine& eng);

    /** 
     * Constructor
     * @param params set of parameters (see parameters.hpp)
     */
    PosteriorModel(size_t dim, Parameters params, randEngine& eng);

    /** 
     * Default destructor
     */
    virtual ~PosteriorModel();

   
    virtual void updateHyperParameters() = 0;
    virtual void fitSurrogateModel() = 0;
    virtual void updateSurrogateModel() = 0;

    virtual double evaluateCriteria(const vectord& query) = 0;
    virtual void updateCriteria(const vectord& query) = 0;

    virtual bool criteriaRequiresComparison() = 0;
    virtual void setFirstCriterium() = 0;
    virtual bool setNextCriterium(const vectord& prevResult) = 0;
    virtual std::string getBestCriteria(vectord& best) = 0;


    void setSamples(const matrixd &x, const vectord &y);
    void setSamples(const matrixd &x);
    void setSamples(const vectord &y);
    void setSample(const vectord &x, double y);
    void addSample(const vectord &x, double y);
    double getValueAtMinimum();
    vectord getPointAtMinimum();

    void plotDataset(TLogLevel level);

    const Dataset* getData();
    virtual ProbabilityDistribution* getPrediction(const vectord& query) = 0;


  protected:
    Parameters mParameters;                     ///< Configuration parameters
    size_t mDims;                                    ///< Number of dimensions
    Dataset mData;                ///< Dataset (x-> inputs, y-> labels/output)
    MeanModel mMean;

  private:
    PosteriorModel();
  };

  inline vectord PosteriorModel::getPointAtMinimum() 
  { return mData.getPointAtMinimum(); };
  
  inline double PosteriorModel::getValueAtMinimum()
  { return mData.getValueAtMinimum(); };

  inline void PosteriorModel::plotDataset(TLogLevel level)
  { mData.plotData(level); }

  inline const Dataset* PosteriorModel::getData()
  { return &mData; }

} //namespace bayesopt


#endif
