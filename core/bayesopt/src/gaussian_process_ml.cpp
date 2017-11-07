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
#include "gaussian_process_ml.hpp"

namespace bayesopt
{

  namespace ublas = boost::numeric::ublas;

  GaussianProcessML::GaussianProcessML(size_t dim, Parameters params, 
				       const Dataset& data, 
				       MeanModel& mean, randEngine& eng):
    HierarchicalGaussianProcess(dim, params, data, mean, eng)
  {
    d_ = new GaussianDistribution(eng);
  }  // Constructor



  GaussianProcessML::~GaussianProcessML()
  {
    delete d_;
  } // Default destructor




  double GaussianProcessML::negativeLogLikelihood()
  {
    /* In this case, they are equivalent */
    return negativeTotalLogLikelihood();
  }


  ProbabilityDistribution* GaussianProcessML::prediction( const vectord &query )
  {
    double kq = computeSelfCorrelation(query);
    vectord kn = computeCrossCorrelation(query);
    vectord phi = mMean.getFeatures(query);
  
    vectord v(kn);
    inplace_solve(mL,v,ublas::lower_tag());

    vectord rq = phi - prod(v,mKF);

    vectord rho(rq);
    inplace_solve(mL2,rho,ublas::lower_tag());
    
    double yPred = inner_prod(phi,mWML) + inner_prod(v,mAlphaF);
    double sPred = sqrt( mSigma * (kq - inner_prod(v,v) 
				   + inner_prod(rho,rho)));

    d_->setMeanAndStd(yPred,sPred);
    return d_;
  }

  void GaussianProcessML::precomputePrediction()
  {
    size_t n = mData.getNSamples();
    size_t p = mMean.nFeatures();

    mKF = trans(mMean.mFeatM);
    inplace_solve(mL,mKF,ublas::lower_tag());

    matrixd FKF = prod(trans(mKF),mKF);
    mL2 = FKF;
    utils::cholesky_decompose(FKF,mL2);

    vectord Ky(mData.mY);
    inplace_solve(mL,Ky,ublas::lower_tag());

    mWML = prod(Ky,mKF);
    utils::cholesky_solve(mL2,mWML,ublas::lower());

    mAlphaF = mData.mY - prod(mWML,mMean.mFeatM);
    inplace_solve(mL,mAlphaF,ublas::lower_tag());
    mSigma = inner_prod(mAlphaF,mAlphaF)/(n-p);
  }

} //namespace bayesopt
