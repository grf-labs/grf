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

#include <cstdlib>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include "ublas_trace.hpp"
#include "ublas_elementwise.hpp"
#include "gauss_distribution.hpp"
#include "gaussian_process_normal.hpp"

namespace bayesopt
{

  namespace ublas = boost::numeric::ublas; 
  
  GaussianProcessNormal::GaussianProcessNormal(size_t dim, 
					       Parameters params, 
					       const Dataset& data, 
					       MeanModel& mean,
					       randEngine& eng):
    HierarchicalGaussianProcess(dim,params,data,mean,eng),
    mW0(params.mean.coef_mean.size()), mInvVarW(params.mean.coef_mean.size()), 
    mD(params.mean.coef_mean.size(),params.mean.coef_mean.size())
  {  
    mSigma = params.sigma_s;
    mW0 = params.mean.coef_mean;
    for (size_t ii = 0; ii < params.mean.coef_std.size(); ++ii)
      {
	double varii = params.mean.coef_std[ii] * params.mean.coef_std[ii];
	mInvVarW(ii) = 1/varii;
      }
     d_ = new GaussianDistribution(eng);
  }  // Constructor



  GaussianProcessNormal::~GaussianProcessNormal()
  {
    delete d_;
  } // Default destructor


  ProbabilityDistribution* 
  GaussianProcessNormal::prediction(const vectord &query)
  {
    const double kq = computeSelfCorrelation(query);
    const vectord phi = mMean.getFeatures(query);

    vectord v = computeCrossCorrelation(query);

    inplace_solve(mL,v,ublas::lower_tag());

    vectord rq = phi - prod(v,mKF);

    vectord rho(rq);
    inplace_solve(mD,rho,ublas::lower_tag());
    
    double yPred = inner_prod(phi,mWMap) + inner_prod(v,mVf);
    double sPred = sqrt( mSigma * (kq - inner_prod(v,v) 
			        + inner_prod(rho,rho)));

    if ((boost::math::isnan(yPred)) || (boost::math::isnan(sPred)))
      {
	throw std::runtime_error("Error in prediction. NaN found.");
      }
					

    d_->setMeanAndStd(yPred,sPred);
    return d_;
  }


  double GaussianProcessNormal::negativeLogLikelihood()
  {
    matrixd KK = computeCorrMatrix();
    const size_t n = KK.size1();
    const size_t p = mMean.nFeatures();
  
    vectord v0 = mData.mY - prod(trans(mMean.mFeatM),mW0);
    matrixd WW = zmatrixd(p,p);  //TODO: diagonal matrix
    utils::add_to_diagonal(WW,mInvVarW);
    matrixd FW = prod(trans(mMean.mFeatM),WW);
    KK += prod(FW,mMean.mFeatM);
    matrixd BB(n,n);
    utils::cholesky_decompose(KK,BB);
    inplace_solve(BB,v0,ublas::lower_tag());
    double zz = inner_prod(v0,v0);

    double lik = 1/(2*mSigma) * zz;
    lik += utils::log_trace(BB);
    return lik;
  }



  void GaussianProcessNormal::precomputePrediction()
  {
    const size_t p = mMean.nFeatures();

    mKF = trans(mMean.mFeatM);
    inplace_solve(mL,mKF,ublas::lower_tag());
    //TODO: make one line
    matrixd DD(p,p);
    DD = prod(trans(mKF),mKF);
    utils::add_to_diagonal(DD,mInvVarW);
    utils::cholesky_decompose(DD,mD);

    vectord vn = mData.mY;
    inplace_solve(mL,vn,ublas::lower_tag());
    mWMap = prod(mMean.mFeatM,vn) + utils::ublas_elementwise_prod(mInvVarW,mW0);
    utils::cholesky_solve(mD,mWMap,ublas::lower());

    mVf = mData.mY - prod(trans(mMean.mFeatM),mWMap);
    inplace_solve(mL,mVf,ublas::lower_tag());

    if (boost::math::isnan(mWMap(0)))
      {
	throw std::runtime_error("Error in precomputed prediction. NaN found.");
      }
  }

} //namespace bayesopt
