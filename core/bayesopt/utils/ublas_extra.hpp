/**  \file ublas_extra.hpp \brief Extra functions for Ublas library */
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

#ifndef __UBLAS_EXTRA_HPP__
#define __UBLAS_EXTRA_HPP__

#include <typeinfo>
#include <boost/numeric/ublas/vector.hpp>

namespace bayesopt
{
  /// Extra utils: math functions, ublas helpers, etc.
  namespace utils
  {
    template<class V, class D>
    void append(V& vect, D element)
    {
      typedef typename V::value_type VD;
      assert(typeid(VD) == typeid(D));
      
      // This method is super inefficient but there seems to be the uBlas style.
      const size_t size = vect.size();
      vect.resize(size+1,true);
      vect(size) = element;
    };

    template<class V, class I>
    void erase(V& vect, I begin)
    {
      typedef typename V::iterator VI;
      assert(typeid(VI) == typeid(I));
     
      for(VI it = begin; it != vect.end()-1; ++it)
	{
	  *it = *(it+1); 
	}
      vect.resize(vect.size()-1);
    };

    template<class M>
    void erase_column(M& mat, size_t pos)
    {
      for(size_t i = pos; i < mat.size2()-1; ++i)
	{
	  column(mat,i) = column(mat,i+1);
	}
      mat.resize(mat.size1(),mat.size2()-1);
    };

    template<class M, class V>
    void add_to_diagonal(M& mat, const V& vec)
    {
      assert(mat.size1()==mat.size2());
      assert(mat.size1()==vec.size());
      const size_t ll = vec.size();
      for(size_t ii = 0; ii < ll; ++ii)
	{
	  mat(ii,ii) += vec(ii);
	}
    };

    boost::numeric::ublas::vector<double> array2vector(const double array[], 
						       const size_t n);

  } //  namespace utils
} //namespace bayesopt

#endif
