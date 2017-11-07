
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


#include "bayesopt/bayesopt.hpp"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
//#include "randgen.hpp"
#include "lhs.hpp"
#include "gridsampling.hpp"
#include "log.hpp"

namespace bayesopt
{
  
  // DiscreteModel::DiscreteModel(const vecOfvec &validSet):
  //   BayesOptBase(), mInputSet(validSet)
  // {} // Constructor


  DiscreteModel::DiscreteModel( const vecOfvec &validSet, 
				Parameters parameters):
    BayesOptBase(validSet[0].size(),parameters), mInputSet(validSet)
  {    
    mDims = mInputSet[0].size();    
  } // Constructor

  DiscreteModel::DiscreteModel(const vectori &categories, 
			       Parameters parameters):
   BayesOptBase(categories.size(),parameters)
  {    
    mDims = categories.size();    
    utils::buildGrid(categories,mInputSet);
  }


  DiscreteModel::~DiscreteModel()
  {} // Default destructor



  // PROTECTED
  vectord DiscreteModel::samplePoint()
  {   
    randInt sample(mEngine, intUniformDist(0,mInputSet.size()-1));
    return mInputSet[sample()];
  };

  void DiscreteModel::findOptimal(vectord &xOpt)
  {
    std::vector<double> critv(mInputSet.size());
    std::transform(mInputSet.begin(),mInputSet.end(),critv.begin(),
		   boost::bind(&DiscreteModel::evaluateCriteria,this,_1));

    xOpt = mInputSet[std::distance(critv.begin(),
			 std::max_element(critv.begin(),critv.end()))];
    
    // xOpt = *mInputSet.begin();
    // double min = evaluateCriteria(xOpt);
    
    // for(vecOfvec::iterator it = mInputSet.begin();
    // 	it != mInputSet.end(); ++it)
    //   {
    // 	double current = evaluateCriteria(*it);
    // 	if (current < min)
    // 	  {
    // 	    xOpt = *it;  
    // 	    min = current;
    // 	  }
    //   }
  }

  //In this case, it is the trivial function
  vectord DiscreteModel::remapPoint(const vectord& x)
  { return x; }

  void DiscreteModel::generateInitialPoints(matrixd& xPoints)
  {

    vecOfvec perms = mInputSet;
    
    // By using random permutations, we guarantee that 
    // the same point is not selected twice
    utils::randomPerms(perms,mEngine);
    
    // vectord xPoint(mInputSet[0].size());
    for(size_t i = 0; i < xPoints.size1(); i++)
    {
        const vectord xP = perms[i];
        row(xPoints,i) = xP;
    }
  }

}  // namespace bayesopt


