/**  \file displaygp.hpp \brief Plots the evolution (nonparametric
     process, criteria or contour plots) of 1D and 2D problems. */
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

#ifndef _DISPLAYGP_HPP_
#define _DISPLAYGP_HPP_

#include "bayesopt/bayesopt.hpp"
#include "matplotpp.h"  

namespace bayesopt
{
  namespace utils
  {      
    enum RunningStatus
      {
	RUN, STEP, STOP, NOT_READY
      };

    class DisplayProblem1D :public MatPlot
    { 
    private:
      RunningStatus status;
      size_t state_ii;
      BayesOptBase* bopt_model;
      std::vector<double> lx,ly;

    public:
      DisplayProblem1D();
      void init(BayesOptBase* bopt, size_t dim);
      void setSTEP();
      void toogleRUN();
      void DISPLAY();
    };

    class DisplayProblem2D :public MatPlot
    { 
    private:
      RunningStatus status;
      size_t state_ii;
      BayesOptBase* bopt_model;
      std::vector<double> lx,ly;
      std::vector<double> cx, cy;
      std::vector<double> solx, soly;
      size_t c_points;
      std::vector<double> cX,cY;
      std::vector<std::vector<double> > cZ;

    public:
      DisplayProblem2D();
      void setSolution(vectord sol);
      void prepareContourPlot();
      void init(BayesOptBase* bopt, size_t dim);
      void setSTEP();
      void toogleRUN();
      void DISPLAY();
    };


  } //namespace utils

} //namespace bayesopt


#endif
