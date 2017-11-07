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
#include "log.hpp"
#include "lhs.hpp"
#include "mcmc_sampler.hpp"

namespace bayesopt
{
  MCMCSampler::MCMCSampler(RBOptimizable* rbo, size_t dim, randEngine& eng):
    mtRandom(eng), obj(new RBOptimizableWrapper(rbo))
  {
    mAlg = SLICE_MCMC;
    mDims = dim;
    nBurnOut = 100;
    nSamples = 10;
    mStepOut = true;
    mSigma = svectord(dim,6);
  };

  MCMCSampler::~MCMCSampler()
  {};

  void MCMCSampler::randomJump(vectord &x)
  {
    randNFloat sample( mtRandom, normalDist(0,1) );
    FILE_LOG(logERROR) << "Doing random jump.";
    for(vectord::iterator it = x.begin(); it != x.end(); ++it)
      {
	*it = sample()*6;
      }
    FILE_LOG(logERROR) << "Likelihood." << x << " | " << obj->evaluate(x);
  }

  //TODO: Include new algorithms when we add them.
  void MCMCSampler::burnOut(vectord &x)
  {
    for(size_t i=0; i<nBurnOut; ++i)  
      {
	try
	  {
	    sliceSample(x);
	  }
	catch(std::runtime_error& e)
	  {
	    FILE_LOG(logERROR) << e.what();
	    randomJump(x);
	  }
      }
  }


  void MCMCSampler::sliceSample(vectord &x)
  {
    randFloat sample( mtRandom, realUniformDist(0,1) );
    size_t n = x.size();

    std::vector<int> perms = utils::return_index_vector(0,n);

    utils::randomPerms(perms, mtRandom);

    for (size_t i = 0; i<n; ++i)
      {
	const size_t ind = perms[i];
	const double sigma = mSigma(ind);

	const double y_max = -obj->evaluate(x);
	const double y = y_max+std::log(sample());  

	if (y == 0.0) 
	  {
	    throw std::runtime_error("Error in MCMC: Initial point"
				     " out of support region."); 
	  }

	// Step out
	const double x_cur = x(ind);
	const double r = sample();
	double xl = x_cur - r * sigma;
	double xr = x_cur + (1-r)*sigma;

	if (mStepOut)
	  {
	    x(ind) = xl;
	    while (-obj->evaluate(x) > y) { x(ind) -= sigma; }
	    xl = x(ind);

	    x(ind) = xr;
	    while (-obj->evaluate(x) > y) { x(ind) += sigma; }
	    xr = x(ind);
	  }

	//Shrink
	bool on_slice = false;
	while (!on_slice)
	  {
	    x(ind) = (xr-xl) * sample() + xl;
	    if (-obj->evaluate(x) < y)
	      {
		if      (x(ind) > x_cur)  xr = x(ind);
		else if (x(ind) < x_cur)  xl = x(ind);
		else throw std::runtime_error("Error in MCMC. Slice colapsed.");
	      }
	    else
	      {
		on_slice = true;
	      }
	  }
      }
  }

  //TODO: Include new algorithms when we add them.
  void MCMCSampler::run(vectord &Xnext)
  {
    randFloat sample( mtRandom, realUniformDist(0,1) );
    if (nBurnOut>0) burnOut(Xnext);

    mParticles.clear();
    for(size_t i=0; i<nSamples; ++i)  
      {
	try
	  {
	    sliceSample(Xnext);
	  }
	catch(std::runtime_error& e)
	  {
	    FILE_LOG(logERROR) << e.what();
	    randomJump(Xnext);
	  }
	mParticles.push_back(Xnext);
      }
    printParticles();
  }

}// namespace bayesopt
