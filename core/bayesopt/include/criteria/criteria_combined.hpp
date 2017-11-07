/**  \file criteria_combined.hpp \brief Abstract module for combined criteria */
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

#ifndef  _CRITERIA_COMBINED_HPP_
#define  _CRITERIA_COMBINED_HPP_

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include "criteria_functors.hpp"


namespace bayesopt
{
  /**\addtogroup CriteriaFunctions */
  //@{

  /// Abstract class for combined criteria functions
  class CombinedCriteria: public Criteria
  {
  public:
    virtual ~CombinedCriteria() {};
    virtual void init(NonParametricProcess *proc) 
    { 
      mProc = proc;
    };

    void pushCriteria(Criteria* crit)
    {
      mCriteriaList.push_back(crit);
    };

    void setParameters(const vectord &theta) 
    {
      using boost::numeric::ublas::subrange;
      const size_t np = mCriteriaList.size();
      vectori sizes(np);

      for (size_t i = 0; i < np; ++i)
	{
	  sizes(i) = mCriteriaList[i].nParameters();
	}

      if (theta.size() != norm_1(sizes))
	{
	  FILE_LOG(logERROR) << "Wrong number of criteria parameters"; 
	  throw std::invalid_argument("Wrong number of criteria parameters");
	}

      size_t start = 0;
      for (size_t i = 0; i < np; ++i)
	{
	  mCriteriaList[i].setParameters(subrange(theta,start,start+sizes(i)));
	  start += sizes(i);
	}
    };

    size_t nParameters() 
    {
      size_t sum = 0;
      for (size_t i = 0; i < mCriteriaList.size(); ++i)
	{
	  sum += mCriteriaList[i].nParameters();
	}
      return sum;
    };

  protected:
    boost::ptr_vector<Criteria> mCriteriaList;
  };

  //@}

} //namespace bayesopt


#endif
