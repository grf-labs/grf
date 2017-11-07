/**  \file trace_ublas.hpp \brief Trace computation for uBlas */
/*
-----------------------------------------------------------------------------
   Copyright (C) 2011 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU Affero General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Affero General Public License for more details.

   You should have received a copy of the GNU Affero General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
-----------------------------------------------------------------------------
*/
#ifndef __TRACE_UBLAS_HPP__
#define __TRACE_UBLAS_HPP__

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace bayesopt
{
  namespace utils
  {

    template<class E>
    typename E::value_type trace(const E &A)
    {
      const size_t n = (std::min)(A.size1(),A.size2());
      typename E::value_type sum = 0;
      for (size_t i=0; i<n; ++i)
	sum += A(i,i);

      return sum; 
    }

    template<class E>
    typename E::value_type log_trace(const E &A)
    {
      const size_t n = (std::min)(A.size1(),A.size2());
      typename E::value_type sum = 0;
      for (size_t i=0; i<n; ++i)
	sum += std::log(A(i,i));

      return sum; 
    }


    template<class E1, class E2>
    typename E1::value_type trace_prod(const E1 & A, const E2 & B )
    {
      namespace ublas = boost::numeric::ublas;

      const size_t n = (std::min)(A.size1(),B.size2());
      typename E1::value_type sum = 0;
      for (size_t i=0; i<n; ++i)
	sum += ublas::inner_prod(ublas::row(A,i),ublas::column(B,i));

      return sum; 
    }
    
  } //  namespace utils
} //namespace bayesopt

#endif
