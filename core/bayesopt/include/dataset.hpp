/** \file dataset.hpp 
    \brief Dataset model */
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


#ifndef __DATASET_HPP__
#define __DATASET_HPP__

#include "log.hpp"
#include "specialtypes.hpp"
#include "ublas_extra.hpp"

namespace bayesopt
{

  /** \addtogroup NonParametricProcesses */
  /**@{*/

  /** \brief Dataset model to deal with the vector (real) based datasets */
  class Dataset
  {
  public:
    Dataset();
    Dataset(const matrixd& x, const vectord& y);
    virtual ~Dataset();

    void setSamples(const matrixd &x, const vectord &y);
    void setSamples(const matrixd &x);
    void setSamples(const vectord &y);
    void addSample(const vectord &x, double y);
    double getSampleY(size_t index) const;
    vectord getSampleX(size_t index) const;
    double getLastSampleY() const;
    vectord getLastSampleX() const;

    void plotData(TLogLevel level);

    vectord getPointAtMinimum() const;
    double getValueAtMinimum() const;
    size_t getNSamples() const;
    void updateMinMax( size_t i );

    vecOfvec mX;                                         ///< Data inputs
    vectord mY;                                          ///< Data values

  private:
    size_t mMinIndex, mMaxIndex;	
  };


  //// Inline methods
  
  inline void Dataset::addSample(const vectord &x, double y)
  {
    mX.push_back(x); utils::append(mY,y);
    updateMinMax(mY.size()-1);
  }

  inline double Dataset::getSampleY(size_t index) const
  { return mY(index);  }

  inline vectord Dataset::getSampleX(size_t index) const
  { return mX[index];  }

  inline double Dataset::getLastSampleY() const
  { return mY(mY.size()-1); }

  inline vectord Dataset::getLastSampleX() const
  { return mX[mX.size()-1]; }


  inline vectord Dataset::getPointAtMinimum() const { return mX[mMinIndex]; };
  inline double Dataset::getValueAtMinimum() const { return mY(mMinIndex); };
  inline size_t Dataset::getNSamples() const { return mY.size(); };
  inline void Dataset::updateMinMax( size_t i )
  {
    if ( mY(mMinIndex) > mY(i) )       mMinIndex = i;
    else if ( mY(mMaxIndex) < mY(i) )  mMaxIndex = i;
  };

  /**@}*/
  
} //namespace bayesopt

#endif
