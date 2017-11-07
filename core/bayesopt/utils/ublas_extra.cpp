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

#include "ublas_extra.hpp"


namespace bayesopt
{
  namespace utils
  {
    boost::numeric::ublas::vector<double> array2vector(const double array[], 
						       const size_t n)
    {
      boost::numeric::ublas::vector<double> v(n);
      std::copy(array, array+n, v.begin());
      return v;
    }

  } //  namespace utils
} //namespace bayesopt
