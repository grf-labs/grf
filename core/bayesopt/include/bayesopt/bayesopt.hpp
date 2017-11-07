/**  \file bayesopt.hpp \brief BayesOpt main C++ interface */
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

#ifndef  _BAYESOPTAPI_HPP_
#define  _BAYESOPTAPI_HPP_

#include "bayesoptbase.hpp"

namespace bayesopt  {

  // Forward declarations
  namespace utils 
  {  
    template <class V> 
    class BoundingBox;
  }
  class NLOPT_Optimization;
  class CritCallback;


  
  /** \addtogroup BayesOpt */
  /**@{*/

  /**
   * \brief Bayesian optimization for functions in continuous input spaces. 
   *
   * It requires box constrains for the input space. More exactly:
   * \f$ f(x):\mathcal{D} \subset \mathbb{R}^n \Rightarrow \mathbb{R} \f$.
   *
   * Usage:
   * \code{.cpp}
   *   ContinuousModel opt(dim,params);
   *   vectord result(dim), lBound(dim), uBound(dim);
   *   \\.. Define bounds
   *   opt.setBoundingBox(lBound,uBound);
   * \endcode  
   *
   * Optimization can be run in batch mode calling
   * \code{.cpp}
   *   opt.optimize(result);
   * \endcode  
   * or step by step.
   * \code{.cpp}
   *   opt.initiliazeOptimization();
   *   \\...
   *   opt.stepOptimization();
   *   \\...
   *   result getFinalResult();
   * \endcode  
   *
   * This model can also be used for discrete/integer data, provided
   * that the callback function provides the corresponding casting or
   * nearest neighbour.
   * 
   * \see BayesOptBase about how to run the optimization
   */
  class BAYESOPT_API ContinuousModel: public BayesOptBase
  {
  public:
    /** 
     * Constructor
     * @param dim number of input dimensions
     * @param params set of parameters (see parameters.h)
     */
    ContinuousModel(size_t dim, Parameters params);

    /**  Default destructor  */
    virtual ~ContinuousModel();

    /** 
     * \brief Sets the bounding box. 
     *
     * @param lowerBound vector with the lower bounds of the hypercube
     * @param upperBound vector with the upper bounds of the hypercube
     */
    void setBoundingBox(const vectord &lowerBound,
			const vectord &upperBound);


  protected:
    /** Sample a single point in the input space. Used for epsilon
	greedy exploration. */
    vectord samplePoint();

    /** 
     * \brief Call the inner optimization method to find the optimal
     * point acording to the criteria.  
     * @param xOpt optimal point
     */
    void findOptimal(vectord &xOpt);

    /** Remap the point x to the original space (e.g.:
	unnormalization) */
    vectord remapPoint(const vectord& x);

    /** Selects the initial set of points to build the surrogate model. */
    void generateInitialPoints(matrixd& xPoints);

  private:
    boost::scoped_ptr<utils::BoundingBox<vectord> > mBB;      ///< Bounding Box (input space limits)
    boost::scoped_ptr<NLOPT_Optimization> cOptimizer;
    boost::scoped_ptr<CritCallback> mCallback;

  private:
    ContinuousModel();                       ///< Default constructor forbidden.
  };
  

  /**
   * \brief Bayesian optimization for functions in discrete spaces. 
   *
   * The discrete space can be created in two ways depending on the constructor used:
   *  -# A set of discrete points (real vectors) in a real space. 
   *  -# A set of categories with different values for each category.
   *
   * The kind of models used in this library are more suitable for
   * problems of the first point. However, it can also be used with
   * problems in the second category. In that case, we recommend to
   * use the Hamming kernel function.
   *
   *
   * Usage:
   * \code{.cpp}
   *   DiscreteModel opt(validSetOfPoints,params);
   *   \\ or
   *   DiscreteModel opt(categories,params);
   * \endcode  
   *
   * Optimization can be run in batch mode calling
   * \code{.cpp}
   *   opt.optimize(result);
   * \endcode  
   * or step by step.
   * \code{.cpp}
   *   opt.initiliazeOptimization();
   *   \\...
   *   opt.stepOptimization();
   *   \\...
   *   result getFinalResult();
   * \endcode  
   *
   * \see BayesOptBase about how to run the optimization
   * \see HammingKernel for categorical data.
   */
  class BAYESOPT_API DiscreteModel : public BayesOptBase
  {
  public:
    /** 
     * Constructor for real-valued discrete data
     * @param validSet  Set of potential inputs
     * @param params set of parameters (see parameters.h)
     */
    DiscreteModel(const vecOfvec &validSet, Parameters params);

    /** 
     * Constructor for categorical data
     * @param number of categories per dimension
     * @param params set of parameters (see parameters.h)
     */
    DiscreteModel(const vectori &categories, Parameters params);
    
    /** Default destructor  */
    virtual ~DiscreteModel();
    
  protected:
    /** Sample a single point in the input space. Used for epsilon
	greedy exploration. */
    vectord samplePoint();

    /** 
     * \brief Call the inner optimization method to find the optimal
     * point acording to the criteria.  
     * @param xOpt optimal point
     */
    void findOptimal(vectord &xOpt);

    /** Remap the point x to the original space  */
    vectord remapPoint(const vectord& x);

    /** Selects the initial set of points to build the surrogate model. */
    void generateInitialPoints(matrixd& xPoints);

  private:
    vecOfvec mInputSet;               ///< List of input points

    DiscreteModel();         ///< Default constructor forbidden.
  };


  /**@}*/


}  //namespace bayesopt


#endif
