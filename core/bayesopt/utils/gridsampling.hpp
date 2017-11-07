/**  \file gridsampling.hpp \brief Regular grid sampling. */
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

#ifndef _GRIDSAMPLING_HPP_
#define _GRIDSAMPLING_HPP_

#include "specialtypes.hpp"

namespace bayesopt
{
  namespace utils
  {      
    
    void deepenGrid(size_t comp, const vectori ndivs, 
		vectord& x, vecOfvec& result)
    {
      if (comp == x.size())
	{
	  result.push_back(x);
	}
      else
	{
	  for (size_t i = 0; i<ndivs(comp); ++i)
	    {
	      x(comp) = static_cast<double>(i);
	      deepenGrid(comp+1,ndivs, x, result);
	    }
	}
    };

    void buildGrid(const vectori& dims, vecOfvec& result)
    {
      if (result.size()){  result.clear(); }

      vectord x(dims.size());
      deepenGrid(0,dims,x,result);
    };


  } //namespace utils

} //namespace bayesopt

#endif
