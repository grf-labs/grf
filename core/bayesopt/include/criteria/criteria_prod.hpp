/**  \file criteria_prod.hpp \brief Product of multiple criteria */
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

#ifndef  _CRITERIA_PROD_HPP_
#define  _CRITERIA_PROD_HPP_

#include "criteria_combined.hpp"


namespace bayesopt
{
  /**\addtogroup CriteriaFunctions */
  //@{

  /// Product of criterion functions. Unintiutive, but it might come
  /// handy in the future.
  class ProdCriteria: public CombinedCriteria
  {
  public:
    virtual ~ProdCriteria() {};
  
    double operator() (const vectord &x)   
    {
      double prod = 1.0;
      for(size_t i = 0; i<mCriteriaList.size(); ++i)
	{ 
	  prod *= mCriteriaList[i](x); 
	}
      return prod;
    };

    std::string name() {return "cProd";};
  };

  //@}

} //namespace bayesopt


#endif
