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
#include <stdexcept>
#include <sstream>
#include <iostream>
#include "parser.hpp"

namespace bayesopt 
{  
  namespace utils 
  {
    /**
     * Parse expresions of the form Parent(Child1, Child2). The "childs"
     * can also be expressions of the same type.
     */
    void parseExpresion(std::string input, std::string& parent,
		       std::string& child1, std::string& child2)
    {
      std::stringstream is(input);
      std::stringstream os(std::stringstream::out);
      std::stringstream os1(std::stringstream::out);
      std::stringstream os2(std::stringstream::out);
      char c;
      int i = 0, j = 0;
      while (is >> c) 
	{
	  if (i < 0) throw std::runtime_error("Error parsing expression:" + input);
	  
	  if (c == ' ') /* pass */;
	  else if (c == '(') i++;
	  else if (c == ')') i--;
	  else if ((i == 1) && (c == ',')) j++;
	  else 
	    {
	      if (i == 0) os << c;
	      else if (j == 0) os1 << c;
	      else os2 << c;
	    }
	}
      if (i != 0) throw std::runtime_error("Error parsing expression:" + input);

      parent = os.str();
      child1 = os1.str();
      child2 = os2.str();
    }

    /**
     * Parse expresions of the form Parent(Child1, ... ,ChildN). The "childs"
     * can also be expressions of the same type.
     */
    void parseExpresion(std::string input, std::string& parent,
		       std::vector<std::string>& childs)
    {
      std::stringstream is(input);
      std::stringstream os(std::stringstream::out);
      std::stringstream os1(std::stringstream::out);
      char c;
      int i = 0;

      childs.clear();
      while (is >> c) 
	{
	  if (i < 0) throw std::runtime_error("Error parsing expression:" + input);

	  if (c == ' ') /* pass */;
	  else if (c == '(') 
	    {
	      ++i; 
	      if (i > 1) os1 << c;
	    }
	  else if ((c == ')') && (--i == 0))
	    {
	      childs.push_back(os1.str());	      
	    }
	  else if ((i == 1) && (c == ',')) 
	    {
	      childs.push_back(os1.str());
	      os1.str(std::string());
	    }
	  else 
	    {
	      if (i == 0) os << c;
	      else os1 << c;
	    }
	}
      if (i != 0) throw std::runtime_error("Error parsing expression:" + input);

      parent = os.str();
    } 
    
    /**
     * Splits the input string with a delimiter to extract elements
     */
    void split(std::string &input, char delim, std::vector<std::string> &elems){
        std::stringstream ss(input);
        std::string item;
        while (std::getline(ss, item, delim)){
            elems.push_back(item);
        }
    }
  } //namespace utils

} //namespace bayesopt

