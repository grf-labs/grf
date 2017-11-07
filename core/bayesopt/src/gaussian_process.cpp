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

#include "ublas_trace.hpp"
#include "gaussian_process.hpp"

namespace bayesopt
{

  namespace ublas = boost::numeric::ublas;

  GaussianProcess::GaussianProcess(size_t dim, Parameters params, 
				   const Dataset& data, MeanModel& mean,
				   randEngine& eng):
    ConditionalBayesProcess(dim, params, data, mean, eng)
  {
    mSigma = params.sigma_s;
    d_ = new GaussianDistribution(eng);
  }  // Constructor


  GaussianProcess::~GaussianProcess()
  {
    delete d_;
  } // Default destructor


  double GaussianProcess::negativeTotalLogLikelihood()
  {
    // In this case it is equivalent.
    return negativeLogLikelihood();
  }


  double GaussianProcess::negativeLogLikelihood()
  {
    const matrixd K = computeCorrMatrix();
    const size_t n = K.size1();
    matrixd L(n,n);
    utils::cholesky_decompose(K,L);

    vectord alpha(mData.mY-mMean.muTimesFeat());
    inplace_solve(L,alpha,ublas::lower_tag());
    double loglik = ublas::inner_prod(alpha,alpha)/(2*mSigma);
    loglik += utils::log_trace(L);
    return loglik;
  }

  ProbabilityDistribution* GaussianProcess::prediction(const vectord &query)
  {
    const double kq = computeSelfCorrelation(query);
    const vectord kn = computeCrossCorrelation(query);
    

    vectord vd(kn);
    inplace_solve(mL,vd,ublas::lower_tag());
    double basisPred = mMean.muTimesFeat(query);
    double yPred = basisPred + ublas::inner_prod(vd,mAlphaV);
    double sPred = sqrt(mSigma*(kq - ublas::inner_prod(vd,vd)));
    
    d_->setMeanAndStd(yPred,sPred);
    return d_;
  }


  void GaussianProcess::precomputePrediction()
  {
    const size_t n = mData.getNSamples();
  
    mAlphaV.resize(n,false);
    mAlphaV = mData.mY-mMean.muTimesFeat();
    inplace_solve(mL,mAlphaV,ublas::lower_tag());
  }
	
} //namespace bayesopt
