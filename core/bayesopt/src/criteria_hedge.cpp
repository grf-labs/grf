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
#include <numeric>
#include <algorithm>
#include "boost/bind.hpp"
#include "log.hpp"
#include "criteria/criteria_hedge.hpp"

namespace bayesopt
{

  /** 
   * \brief Softmax function
   * 
   * @param g gain function
   * @param eta smoothness coefficient
   * @return 
   */
  inline double softmax(double g, double eta) 
  {
    return exp(eta*g);
  };

  //////////////////////////////////////////////////////////////////////
  GP_Hedge::GP_Hedge(){};

  void GP_Hedge::init(NonParametricProcess *proc) 
  { 
    mProc = proc;

    size_t n = mCriteriaList.size();
    if (!n)
      {
	throw std::logic_error("Criteria list should be created (pushed)"
			       " before initializing combined criterion.");
      }

    loss_ = zvectord(n); 
    gain_ = zvectord(n); 
    prob_ = zvectord(n);
    cumprob_ = zvectord(n);
  };

  void GP_Hedge::initialCriteria()
  {
    mIndex = 0;
    mCurrentCriterium = &mCriteriaList[mIndex];
    mBestLists.clear();
  };

  bool GP_Hedge::rotateCriteria()
  { 
    ++mIndex;
    if (mIndex >= mCriteriaList.size())
      {
	return false;
      }
    else
      {
	mCurrentCriterium = &mCriteriaList[mIndex];
	return true;
      }
  };
  
  void GP_Hedge::pushResult(const vectord& prevResult)
  {
    loss_(mIndex) = computeLoss(prevResult);
    mBestLists.push_back(prevResult);    
  };
  
  std::string GP_Hedge::getBestCriteria(vectord& best)
  { 
    int optIndex = update_hedge();
    best = mBestLists[optIndex];
    return mCriteriaList[optIndex].name();
  };


  int GP_Hedge::update_hedge()
  {
    // We just care about the differences
    double max_l = *std::max_element(loss_.begin(),loss_.end());
    loss_ += svectord(loss_.size(),max_l);

    // To avoid overflow
    double mean_g = std::accumulate(gain_.begin(),gain_.end(),0.0) 
      / static_cast<double>(gain_.size());
    gain_ -= svectord(gain_.size(),mean_g);

    // Optimal eta according to Shapire
    double max_g = *std::max_element(gain_.begin(),gain_.end());
    double eta = (std::min)(10.0,sqrt(2.0*log(3.0)/max_g));
    
    // Compute probabilities
    std::transform(gain_.begin(), gain_.end(), prob_.begin(),
		   boost::bind(softmax,_1,eta));
       
    //Normalize
    double sum_p =std::accumulate(prob_.begin(),prob_.end(),0.0);
    prob_ /= sum_p;

    //Update bandits gain
    gain_ -= loss_;

    std::partial_sum(prob_.begin(), prob_.end(), cumprob_.begin(), 
		     std::plus<double>());

    randFloat sampleUniform( *mtRandom, realUniformDist(0,1));
    double u = sampleUniform();

    for (size_t i=0; i < cumprob_.size(); ++i)
      {
	if (u < cumprob_(i))
	  return i;
      }
    FILE_LOG(logERROR) << "Error updating Hedge algorithm. " 
		       << "Selecting first criteria by default.";
    return 0;
  };


} //namespace bayesopt
