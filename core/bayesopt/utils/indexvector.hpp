/**  \file indexvector.hpp \brief Generators for index vectors. */
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

#ifndef __INDEX_VECTOR_HPP__
#define __INDEX_VECTOR_HPP__

#include <algorithm>

namespace bayesopt
{
  namespace utils
  {
    /**
     * \brief Simple class to generate sequences of unique numbers
     */
    class CUnique 
    {
    private:
      int current;
    public:
      CUnique(int initial = 1) {current=initial-1;}
      int operator()() {return ++current;}
    };

    /** 
     * Generates a vector of indexes (1..n)
     * @param n vector size
     * @return index vector
     */
    inline std::vector<int> return_index_vector(size_t n)
    {
      CUnique UniqueNumber;
      std::vector<int> arr(n);
      generate (arr.begin(), arr.end(), UniqueNumber);
      return arr;
    };

    /** 
     * Generates a vector of indexes starting at A and size_t n
     * @param a starting point
     * @param n vector size
     * @return index vector
     */
    inline std::vector<int> return_index_vector(int a, size_t n)
    {
      CUnique UniqueNumber(a);
      std::vector<int> arr(n);
      generate (arr.begin(), arr.end(), UniqueNumber);
      return arr;
    };


    /** 
     * Modify a vector of indexes (0..n)
     * @param arr vector
     */
    inline void modify_index_vector(std::vector<int>& arr)
    {
      CUnique UniqueNumber;
      generate (arr.begin(), arr.end(), UniqueNumber);
    };

  } //namespace utils
} //namespace bayesopt

#endif
