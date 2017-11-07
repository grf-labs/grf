
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


#include "conditionalbayesprocess.hpp"
#include "log.hpp"
//#include "optimizekernel.hpp"	


namespace bayesopt
{
  ConditionalBayesProcess::ConditionalBayesProcess(size_t dim, Parameters parameters, 
						   const Dataset& data, 
						   MeanModel& mean, randEngine& eng):
    KernelRegressor(dim,parameters,data,mean,eng)
  {}

  ConditionalBayesProcess::~ConditionalBayesProcess()
  {}


  double ConditionalBayesProcess::evaluateKernelParams()
  { 
    switch(mScoreType)
      {
      case SC_MTL:
	return negativeTotalLogLikelihood();
      case SC_ML:
	return negativeLogLikelihood();
      case SC_MAP:
	// It is a minus because the prior is the positive and we want
	// the negative.
	return negativeLogLikelihood()-mKernel.kernelLogPrior();
      case SC_LOOCV:
	return negativeCrossValidation(); 
      default:
	throw std::invalid_argument("Learning type not supported");
      }	  
  }


  double ConditionalBayesProcess::negativeCrossValidation()
  {
    // This is highly ineffient implementation for comparison purposes.
    Dataset data(mData);

    size_t n = data.getNSamples();
    double sum = 0.0;

    matrixd tempF(mMean.mFeatM);


    // We take the first element, use it for validation and then paste
    // it at the end. Thus, after every iteration, the first element
    // is different and, at the end, all the elements should have
    // rotated.
    for(size_t i = 0; i<n; ++i)
      {
	// Take the first element
	const double y = data.getSampleY(0);
	const vectord x = data.getSampleX(0);

	// Remove it for cross validation
	data.mX.erase(data.mX.begin()); 
	utils::erase(data.mY,data.mY.begin());
	utils::erase_column(mMean.mFeatM,0);

	// Compute the cross validation
	computeCholeskyCorrelation();
	precomputePrediction(); 
	ProbabilityDistribution* pd = prediction(x);
	sum += std::log(pd->pdf(y));

	//Paste it back at the end
	data.addSample(x,y);
	mMean.mFeatM.resize(mMean.mFeatM.size1(),mMean.mFeatM.size2()+1);  
	mMean.mFeatM = tempF;
      }
    std::cout << "End" << data.getNSamples();
    return -sum;   //Because we are minimizing.
  }

} // namespace bayesopt
