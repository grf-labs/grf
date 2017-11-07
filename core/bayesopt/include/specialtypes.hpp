/**  \file specialtypes.hpp \brief Boost vector and matrix types */
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

#ifndef __SPECIALTYPES_HPP__
#define __SPECIALTYPES_HPP__

#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

typedef boost::numeric::ublas::vector<int>                      vectori;
typedef boost::numeric::ublas::vector<double>                   vectord;
typedef boost::numeric::ublas::zero_vector<double>             zvectord;
typedef boost::numeric::ublas::scalar_vector<double>           svectord;
typedef boost::numeric::ublas::matrix<double>                   matrixd;
typedef boost::numeric::ublas::zero_matrix<double>             zmatrixd;

typedef std::vector<vectord>                                   vecOfvec;

// Surprisingly, this is the most efficient version of a growing
// matrix for uBlas, but I leave here the old experiments because it
// might change in the future.
typedef boost::numeric::ublas::matrix<double>                 covMatrix;

// typedef boost::numeric::ublas::bounded_matrix<double, MAX_ITERATIONS, 
//                                            MAX_ITERATIONS> covMatrix; 

// typedef boost::numeric::ublas::symmetric_matrix<double,lower,row_major,
//                          bounded_array<double,90000> >       covMatrix;

#endif
