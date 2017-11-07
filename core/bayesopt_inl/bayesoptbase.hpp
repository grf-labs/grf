
/**  \file bayesoptbase.hpp \brief BayesOpt common module for interfaces */
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


#ifndef  _BAYESOPTBASE_HPP_
#define  _BAYESOPTBASE_HPP_

#include <boost/scoped_ptr.hpp>
#include <boost/random.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "bayesopt/parameters.hpp"


/**
 * Namespace of the library interface
 */
namespace bayesopt {


  //Forward declaration
  class PosteriorModel;
  class ProbabilityDistribution;
  class Dataset;
  class BOptState;

  typedef boost::numeric::ublas::vector<double>                   vectord;
  typedef boost::numeric::ublas::vector<int>                      vectori;
  typedef boost::numeric::ublas::matrix<double>                   matrixd;
  typedef std::vector<vectord>                                   vecOfvec;

  /** \addtogroup BayesOpt
   *  \brief Main module for Bayesian optimization
   */
  /*@{*/

  /**
   * \brief Abstract module for Bayesian optimization.
   *
   * This module provides Bayesian optimization using different
   * non-parametric processes (Gaussian process or Student's t
   * process) as distributions over surrogate functions.
   *
   * \see ContinuousModel for implementations of this module for
   * a continuous input spaces
   *
   * \see DiscreteModel for implementations of this module for
   * a discrete input spaces or categorical input variables
   */
  class BAYESOPT_API BayesOptBase
  {
  public:
    /** 
     * Constructor
     * @param params set of parameters (see parameters.hpp)
     */
    BayesOptBase(size_t dim, Parameters params);

    /** 
     * Default destructor
     */
    virtual ~BayesOptBase();

    /** 
     * \brief Function that defines the actual function to be optimized.
     * This function must be modified (overriden) according to the
     * specific problem.
     *
     * @param query point to be evaluated. 
     * @return value of the function at the point evaluated.
     */
    virtual double evaluateSample( const vectord &query ) = 0;
    

    /** 
     * \brief This function checks if the query is valid or not. It can
     * be used to introduce arbitrary constrains. Since the Gaussian
     * process assumes smoothness, constrains are managed by the inner
     * optimizer (e.g.:DIRECT), being highly time consuming. If the
     * constrain is very tricky, DIRECT will need much more function
     * evaluations.
     *
     * Note: This function is experimental. Thus it is not made pure virtual.
     * Using it is completely optional.
     * 
     * @param query point to be evaluated.
     * 
     * @return boolean value showing if the the function is valid at
     *         the query point or not.
     */ 
    virtual bool checkReachability( const vectord &query )
    { return true; };


    /** 
     * \brief Execute the optimization process of the function defined
     * in evaluateSample.
     * 
     * @see evaluateSample
     * @see checkReachability
     *
     * @param bestPoint returns point with the optimum value in a ublas::vector.
     */
    void optimize(vectord &bestPoint);

    /** 
     * \brief Execute ONE step the optimization process of the
     * function defined in evaluateSample.  
     */  
    void stepOptimization();

    /** Initialize the optimization process.  */
    void initializeOptimization();
    
    /** Once the optimization has been perfomed, return the optimal point. */
    vectord getFinalResult();

    /** Saves the current state of the optimization process into a state class. */
    void saveOptimization(BOptState &state);
    
    /** Restores the optimization process of a previous execution */
    void restoreOptimization(BOptState state);

    // Getters and Setters
    ProbabilityDistribution* getPrediction(const vectord& query);
    const Dataset* getData();
    Parameters* getParameters();
    double getValueAtMinimum();
    size_t getCurrentIter();
    double evaluateCriteria(const vectord& query);

  protected:
    /** Get optimal point in the inner space (e.g.: [0-1] hypercube) */
    vectord getPointAtMinimum();

    /** Wrapper for the target function adding any preprocessing or
	constraint. It also maps the box constrains to the [0,1]
	hypercube if applicable. */
    double evaluateSampleInternal( const vectord &query );

    /** Sample a single point in the input space. Used for epsilon
	greedy exploration. */
    virtual vectord samplePoint() = 0;

    /** 
     * \brief Call the inner optimization method to find the optimal
     * point acording to the criteria.  
     * @param xOpt optimal point
     */
    virtual void findOptimal(vectord &xOpt) = 0;
  
    /** Remap the point x to the original space (e.g.:
	unnormalization) */
    virtual vectord remapPoint(const vectord& x) = 0;

    /** Selects the initial set of points to build the surrogate model. */
    virtual void generateInitialPoints(matrixd& xPoints) = 0;

    /** 
     * \brief Print data for every step according to the verbose level
     * 
     * @param iteration iteration number 
     * @param xNext next point
     * @param yNext function value at next point
     */
    void plotStepData(size_t iteration, const vectord& xNext,
			      double yNext);
        
    /** Eases the process of saving a state during initial samples */
    void saveInitialSamples(matrixd xPoints);
    void saveResponse(double yPoint, bool clear);

  protected:
    Parameters mParameters;                    ///< Configuration parameters
    size_t mDims;                                   ///< Number of dimensions
    size_t mCurrentIter;                        ///< Current iteration number
    boost::mt19937 mEngine;                      ///< Random number generator

  private:
    boost::scoped_ptr<PosteriorModel> mModel;
    double mYPrev;
    size_t mCounterStuck;
  private:

    BayesOptBase();

    /** 
     * \brief Selects the next point to evaluate according to a certain
     * criteria or metacriteria
     * 
     * @return next point to evaluate
     */
    vectord nextPoint();  

  };

  /**@}*/



} //namespace bayesopt


#endif
