/*
-----------------------------------------------------------------------------
   Copyright (C) 2011 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU Affero General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Affero General Public License for more details.

   You should have received a copy of the GNU Affero General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
-----------------------------------------------------------------------------
*/

#include "ublas_trace.hpp"
#include "student_t_process_jef.hpp"

namespace bayesopt
{

  namespace ublas = boost::numeric::ublas; 

  StudentTProcessJeffreys::StudentTProcessJeffreys(size_t dim, 
						   Parameters params, 
						   const Dataset& data, 
						   MeanModel& mean,
						   randEngine& eng):
    HierarchicalGaussianProcess(dim, params, data, mean, eng)
  {
    d_ = new StudentTDistribution(eng);
  }  // Constructor



  StudentTProcessJeffreys::~StudentTProcessJeffreys()
  {
    delete d_;
  } // Default destructor




  double StudentTProcessJeffreys::negativeLogLikelihood()
  {
    /* In this case, they are equivalent */
    return negativeTotalLogLikelihood();
  }


  ProbabilityDistribution* 
  StudentTProcessJeffreys::prediction(const vectord &query )
  {
    clock_t start = clock();
    double kq = computeSelfCorrelation(query);
    //    vectord kn = computeCrossCorrelation(query);
    mKernel.computeCrossCorrelation(mData.mX,query,mKn);
    vectord phi = mMean.getFeatures(query);
  
    //    vectord v(mKn);
    inplace_solve(mL,mKn,ublas::lower_tag());

    vectord rho = phi - prod(mKn,mKF);

    //    vectord rho(rq);
    inplace_solve(mL2,rho,ublas::lower_tag());
    
    double yPred = inner_prod(phi,mWML) + inner_prod(mKn,mAlphaF);
    double sPred = sqrt( mSigma * (kq - inner_prod(mKn,mKn) 
				   + inner_prod(rho,rho)));

    d_->setMeanAndStd(yPred,sPred);
    return d_;
  }

  void StudentTProcessJeffreys::precomputePrediction()
  {
    size_t n = mData.getNSamples();
    size_t p = mMean.nFeatures();

    mKn.resize(n);

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
    
    d_->setDof(n-p);  
  }

} //namespace bayesopt
