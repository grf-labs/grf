
/**  \file bopt_state.hpp \brief Representation of a optimization state */
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


#ifndef  _BOPT_STATE_HPP_
#define  _BOPT_STATE_HPP_

#include "bayesopt/parameters.hpp"
#include "specialtypes.hpp"
#include "fileparser.hpp"


/**
 * Namespace of the library interface
 */
namespace bayesopt {

  /**
   * \brief Class that represents the state of an optimization
   * 
   * This class fields stores the state of an optimization, and
   * provides member functions to allow saving the state into a file
   * or reading the state from a file.
   */
    class BOptState{
    public:        
        /** Constructor, parameters set into default */
        BOptState();
        
        /** Creates or overwrites the provided file with the state */
        void saveToFile(std::string filename);
        
        /** Loads the state from the provided file and takes program_params values if needed */
        bool loadFromFile(std::string filename, Parameters &program_params);
        
        // BayesOptBase members
        size_t mCurrentIter;
        size_t mCounterStuck;
        double mYPrev;
        Parameters mParameters;
        
        // Optimization samples
        vecOfvec mX;
        vectord mY;
    private:
        void loadOrSave(utils::FileParser &fp);
    };
} //namespace bayesopt


#endif
