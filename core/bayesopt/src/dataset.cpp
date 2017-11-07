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

#include "dataset.hpp"

#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace bayesopt
{

  Dataset::Dataset(): mMinIndex(0), mMaxIndex(0) {};

  Dataset::Dataset(const matrixd& x, const vectord& y):
    mMinIndex(0), mMaxIndex(0)
  {
    setSamples(x,y);
  };

  Dataset::~Dataset(){};

  void Dataset::setSamples(const matrixd &x, const vectord &y)
  {
    // WARNING: It assumes mX is empty
    mY = y;
    for (size_t i=0; i<x.size1(); ++i)
      {
	mX.push_back(row(x,i));
	updateMinMax(i);
      } 
  };

  void Dataset::setSamples(const vectord &y)
  {
    mY = y;
    for (size_t i=0; i<y.size(); ++i)
      {
	updateMinMax(i);
      } 
  };


  void Dataset::setSamples(const matrixd &x)
  {
    for (size_t i=0; i<x.size1(); ++i)
      {
	mX.push_back(row(x,i));
      } 
  };

  void Dataset::plotData(TLogLevel level)
  {
    // For logging purpose
    FILE_LOG(level) << "Initial points:" ;
    for(size_t i = 0; i < mY.size(); i++)
      {
	FILE_LOG(level) << "X:" << mX[i]
			<< "|Y:" << mY(i);
      }
    const double yPoint = getValueAtMinimum();
    const vectord xPoint = getPointAtMinimum();
    FILE_LOG(level) << "Best point so far:" ;
    FILE_LOG(level) << "X:" << xPoint
		    << "|Y:" << yPoint;
    
  } // plotData


} //namespace bayesopt
