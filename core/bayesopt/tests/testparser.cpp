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

#include <string>
#include <iostream>
#include "parser.hpp"

int main()
{
  std::string test = "One(Two, Three, Four(Five, Six))";
  std::string one;
  std::vector<std::string> vs;
  bayesopt::utils::parseExpresion(test,one,vs);

  std::cout << one << std::endl;
  for(std::vector<std::string>::iterator it = vs.begin();
      it != vs.end(); ++it)
    std::cout << *it << std::endl;
 
  return 0;
}
