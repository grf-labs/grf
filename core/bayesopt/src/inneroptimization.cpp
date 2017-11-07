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
#include <cmath>
#include <nlopt.hpp>
#include "bayesopt/parameters.h"
#include "log.hpp"
#include "inneroptimization.hpp"

namespace bayesopt
{
  void checkNLOPTerror(nlopt_result errortype)
  {
    switch(errortype)
      {
      case -1: FILE_LOG(logERROR) << "NLOPT: General failure"; break;
      case -2: FILE_LOG(logERROR) << "NLOPT: Invalid arguments. Check bounds."; break;
      case -3: FILE_LOG(logERROR) << "NLOPT: Out of memory"; break;
      case -4: FILE_LOG(logERROR) << "NLOPT Warning: Potential roundoff error. " 
				  << "In general, this can be ignored."; break;
      case -5: FILE_LOG(logERROR) << "NLOPT: Force stop."; break;
      default: ;
      }
  }

  const size_t MAX_INNER_EVALUATIONS = 500;   /**< Used per dimmension */

  typedef double (*eval_func)(unsigned int n, const double *x,
			      double *gradient, /* NULL if not needed */
			      void *func_data);


  double run_nlopt(nlopt::algorithm algo, eval_func fpointer,
		   vectord& Xnext, int maxf, const std::vector<double>& vd, 
		   const std::vector<double>& vu, void* objPointer)
  {
    double fmin = 0.0;
    size_t n = Xnext.size(); 
    nlopt::opt opt (algo,n);

    std::vector<double> xstd(n);
    opt.set_lower_bounds(vd);
    opt.set_upper_bounds(vu);
    opt.set_min_objective(fpointer, objPointer);
    opt.set_maxeval(maxf);
    
    // It seems BOBYQA can be unstable if the same point is repeated
    // tested over and over. NLOPT bug?
    opt.set_ftol_rel(1e-12);	
    opt.set_ftol_abs(1e-12);

    std::copy(Xnext.begin(),Xnext.end(),xstd.begin());
      
    try 
      { 
	opt.optimize(xstd, fmin);  
      }
    catch (nlopt::roundoff_limited& e)
      {
	FILE_LOG(logDEBUG) << "NLOPT Warning: Potential roundoff error. " 
			   << "In general, this can be ignored.";
      }

    std::copy(xstd.begin(),xstd.end(),Xnext.begin());
    return fmin;
  }


  NLOPT_Optimization::NLOPT_Optimization(RBOptimizable* rbo, size_t dim):
  mDown(dim),mUp(dim)
  { 
    rbobj = new RBOptimizableWrapper(rbo);       rgbobj = NULL;
    alg = DIRECT;                  maxEvals = MAX_INNER_EVALUATIONS;
    setLimits(zvectord(dim),svectord(dim,1.0));  
  };

  NLOPT_Optimization::NLOPT_Optimization(RGBOptimizable* rgbo, size_t dim):
  mDown(dim),mUp(dim)
  { 
    rbobj = NULL;             rgbobj = new RGBOptimizableWrapper(rgbo);
    alg = DIRECT;             maxEvals = MAX_INNER_EVALUATIONS;
    setLimits(zvectord(dim),svectord(dim,1.0));  
  };

  NLOPT_Optimization::~NLOPT_Optimization()
  {
    if(rbobj != NULL) delete rbobj;
    if(rgbobj != NULL) delete rgbobj;
  }

  double NLOPT_Optimization::localTrialAround(vectord& Xnext)
  {
    assert(mDown.size() == Xnext.size());
    assert(mUp.size() == Xnext.size());
    const size_t n = Xnext.size();

    for (size_t i = 0; i < n; ++i) 
      {
	if (Xnext(i) < mDown[i] || Xnext(i) > mUp[i])
	  {
	    FILE_LOG(logDEBUG) << Xnext;
	    throw std::invalid_argument("Local trial withour proper"
					" initial point.");
	  }
      }

    nlopt::algorithm algo = nlopt::LN_BOBYQA;
    eval_func fpointer = &(NLOPT_Optimization::evaluate_nlopt);
    void* objPointer = static_cast<void *>(rbobj);
    const size_t nIter = 20;
    // std::vector<double> vd(n);
    // std::vector<double> vu(n);

    // for (size_t i = 0; i < n; ++i) 
    //   {
    // 	vd[i] = Xnext(i) - 0.01;
    // 	vu[i] = Xnext(i) + 0.01;
    //   }

    vectord start = Xnext;

    double fmin = run_nlopt(algo,fpointer,Xnext,nIter,
			    mDown,mUp,objPointer);

    FILE_LOG(logDEBUG) << "Near trial " << nIter << "|" 
		       << start << "-> " << Xnext << " f() ->" << fmin;
    
    return fmin;

  }

  double NLOPT_Optimization::run(vectord &Xnext)
  {   
    assert(mDown.size() == Xnext.size());
    assert(mUp.size() == Xnext.size());

    eval_func fpointer;
    void *objPointer;

    size_t n = Xnext.size();
    double fmin = 1;
    int maxf1 = maxEvals*n;
    int maxf2 = 0;    // For a second pass
    const double coef_local = 0.1;
    //int ierror;

    // If Xnext is outside the bounding box, maybe it is undefined
    for (size_t i = 0; i < n; ++i) 
      {
	if (Xnext(i) < mDown[i] || Xnext(i) > mUp[i])
	  {
	    Xnext(i)=(mDown[i]+mUp[i])/2.0;
	  }
      }

    //    nlopt_opt opt;
    nlopt::algorithm algo;
    switch(alg)
      {
      case DIRECT: // Pure global. No gradient
	algo = nlopt::GN_DIRECT_L;
	fpointer = &(NLOPT_Optimization::evaluate_nlopt);
	objPointer = static_cast<void *>(rbobj);
	break;
      case COMBINED: // Combined local-global (80% DIRECT -> 20% BOBYQA). No gradient
	algo = nlopt::GN_DIRECT_L;
	maxf2 = static_cast<int>(static_cast<double>(maxf1)*coef_local);
	maxf1 -= maxf2;  // That way, the number of evaluations is the same in all methods.
	fpointer = &(NLOPT_Optimization::evaluate_nlopt);
	objPointer = static_cast<void *>(rbobj);
	break;
      case BOBYQA:  // Pure local. No gradient
	algo = nlopt::LN_BOBYQA;
	fpointer = &(NLOPT_Optimization::evaluate_nlopt);
	objPointer = static_cast<void *>(rbobj);
	break;
      case LBFGS:  // Pure local. Gradient based
	algo = nlopt::LD_LBFGS;
	fpointer = &(NLOPT_Optimization::evaluate_nlopt_grad);
	objPointer = static_cast<void *>(rgbobj);
	break;
      default: 
	throw std::invalid_argument("Inner optimization algorithm"
				    " not supported");
      }

    if (objPointer == NULL)
      {
	throw std::invalid_argument("Wrong object model "
				    "(gradient/no gradient)");
      }

    fmin = run_nlopt(algo,fpointer,Xnext,maxf1,
		     mDown,mUp,objPointer);

    FILE_LOG(logDEBUG) << "1st opt " << maxf1 << "-> " << Xnext 
		       << " f() ->" << fmin;
    if (maxf2)
      {
	//If the point is exactly at the limit, we may have trouble.
    	for (size_t i = 0; i < n; ++i) 
	  {
	    if (Xnext(i)-mDown[i] < 0.0001)
	      {
		Xnext(i) += 0.0001;
		FILE_LOG(logDEBUG) << "Hacking point for BOBYQA. THIS SHOULD NOT HAPPEN";
	      }
	    if (mUp[i] - Xnext(i) < 0.0001)
	      {
		Xnext(i) -= 0.0001;
		FILE_LOG(logDEBUG) << "Hacking point for BOBYQA. THIS SHOULD NOT HAPPEN";
	      }
	  }

	// BOBYQA may fail in this point. Could it be that EI is not twice differentiable?
	// fmin = run_nlopt(nlopt::LN_BOBYQA,fpointer,Xnext,maxf2,
	// 		 mDown,mUp,objPointer);
	fmin = run_nlopt(nlopt::LN_COBYLA,fpointer,Xnext,maxf2,
			 mDown,mUp,objPointer);
	FILE_LOG(logDEBUG) << "2nd opt " << maxf2 << "-> " << Xnext 
			   << " f() ->" << fmin;
      }

    return fmin;

  } // innerOptimize (uBlas)

  double NLOPT_Optimization::evaluate_nlopt (unsigned int n, const double *x,
					     double *grad, void *my_func_data)

  {
    vectord vx(n);
    std::copy(x,x+n,vx.begin());

    void *objPointer = my_func_data;
    RBOptimizableWrapper* OPTIMIZER = static_cast<RBOptimizableWrapper*>(objPointer);

    return OPTIMIZER->evaluate(vx);
  } /* evaluate_criteria_nlopt */

  double NLOPT_Optimization::evaluate_nlopt_grad (unsigned int n, const double *x,
						  double *grad, void *my_func_data)

  {
    vectord vx(n);
    std::copy(x,x+n,vx.begin());
    
    void *objPointer = my_func_data;
    RGBOptimizableWrapper* OPTIMIZER = static_cast<RGBOptimizableWrapper*>(objPointer);
    

    vectord vgrad = zvectord(n);
    double f =  OPTIMIZER->evaluate(vx,vgrad);
    if (grad && n)  std::copy(vgrad.begin(),vgrad.end(),grad);

    return f;
  } /* evaluate_criteria_nlopt */



}// namespace bayesopt

