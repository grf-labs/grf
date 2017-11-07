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


#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <boost/lexical_cast.hpp>

#include "kernelregressor.hpp"

#include "log.hpp"
#include "ublas_extra.hpp"


namespace bayesopt
{
  KernelRegressor::KernelRegressor(size_t dim, Parameters parameters,
				   const Dataset& data,
				   MeanModel& mean, randEngine& eng):
    NonParametricProcess(dim,parameters,data,mean,eng), 
    mRegularizer(parameters.noise),
    mKernel(dim, parameters),
    mScoreType(parameters.sc_type),
    mLearnType(parameters.l_type),
    mLearnAll(parameters.l_all)
  { }

  KernelRegressor::~KernelRegressor(){}



  void KernelRegressor::updateSurrogateModel()
  {
    const vectord lastX = mData.getLastSampleX();
    vectord newK = computeCrossCorrelation(lastX);
    newK(newK.size()-1) += mRegularizer;   // We add it to the last element
    utils::cholesky_add_row(mL,newK);
    precomputePrediction(); 
  } // updateSurrogateModel


  void KernelRegressor::computeCholeskyCorrelation()
  {
    size_t nSamples = mData.getNSamples();
    mL.resize(nSamples,nSamples);
  
    //  const matrixd K = computeCorrMatrix();
    matrixd K(nSamples,nSamples);
    computeCorrMatrix(K);
    size_t line_error = utils::cholesky_decompose(K,mL);
    if (line_error) 
      {
	throw std::runtime_error("Cholesky decomposition error at line " + 
				 boost::lexical_cast<std::string>(line_error));
      }
  }

  matrixd KernelRegressor::computeDerivativeCorrMatrix(int dth_index)
  {
    const size_t nSamples = mData.getNSamples();
    matrixd corrMatrix(nSamples,nSamples);
    mKernel.computeDerivativeCorrMatrix(mData.mX,corrMatrix,dth_index);
    return corrMatrix;
  }

} //namespace bayesopt
