/** -*- c++ -*- \file bayesoptbase.hpp \brief Bayesian optimization module */
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


#ifndef  _DLL_STUFF_HPP_
#define  _DLL_STUFF_HPP_


/* WINDOWS DLLs stuff */
#if defined (BAYESOPT_DLL) && (defined(_WIN32) || defined(__WIN32__)) && !defined(__LCC__)
  #if defined(bayesopt_EXPORTS)
    #define  BAYESOPT_API __declspec(dllexport)
  #else
    #define  BAYESOPT_API __declspec(dllimport)
  #endif /* MyLibrary_EXPORTS */
#else /* defined (_WIN32) */
 #define BAYESOPT_API
#endif

        
        
#endif
