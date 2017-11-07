
/**  \file param_loader.hpp \brief Allows to load parameters from file */
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


#ifndef  _PARAM_LOADER_HPP_
#define  _PARAM_LOADER_HPP_

#include <fstream>
#include <iostream>
#include <string.h>

#include "bayesopt/parameters.hpp"
#include "specialtypes.hpp"
#include "fileparser.hpp"


/**
 * Namespace of the library interface
 */
namespace bayesopt {
    namespace utils{
        class ParamLoader{
        public:
            /* Loads params from provided file, returns true if file exists */
            static bool load(std::string filename, Parameters &par);
            /* This one could be useful to create param files from existing hard-coded params */
            static void save(std::string filename, Parameters &par); 
        private:
            ParamLoader(){}
        
            /* private function to be used by load and save functions */
            static void loadOrSave(utils::FileParser &fp, Parameters &par);   
        };
    }
} //namespace bayesopt


#endif
