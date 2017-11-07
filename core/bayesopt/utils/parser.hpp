/**  \file parser.hpp \brief Functions to parse strings */
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

#ifndef  _PARSER_HPP_
#define  _PARSER_HPP_

#include <vector>

namespace bayesopt 
{  
  namespace utils 
  {
    /**
     * Parse expresions of the form Parent(Child1, Child2). The "childs"
     * can also be expressions of the same type.
     */
    void parseExpresion(std::string input, std::string& parent,
			std::string& child1, std::string& child2);

    /**
     * Parse expresions of the form Parent(Child1, ... ,ChildN). The "childs"
     * can also be expressions of the same type.
     */
    void parseExpresion(std::string input, std::string& parent,
		       std::vector<std::string>& childs);

    /**
     * Splits the input string with a delimiter to extract elements
     */
    void split(std::string &input, char delim, std::vector<std::string> &elems);
    
  } //namespace utils

} //namespace bayesopt

#endif
