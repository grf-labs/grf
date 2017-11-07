/**  \file randgen.hpp 
    \brief Boost types for random number generation */
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

#ifndef  _RANDGEN_HPP_
#define  _RANDGEN_HPP_

//#include <boost/version.hpp>
#include <boost/random.hpp>

// Types for pseudorandom number generators.

typedef boost::mt19937                                              randEngine;

typedef boost::uniform_real<>				       realUniformDist;
typedef boost::uniform_int<> 					intUniformDist;
typedef boost::normal_distribution<>                                normalDist;
typedef boost::gamma_distribution<>                                  gammaDist;

typedef boost::variate_generator<randEngine&, intUniformDist>          randInt;
typedef boost::variate_generator<randEngine&, normalDist>           randNFloat;
typedef boost::variate_generator<randEngine&, gammaDist>            randGFloat;
typedef boost::variate_generator<randEngine&, realUniformDist>       randFloat;

#endif
