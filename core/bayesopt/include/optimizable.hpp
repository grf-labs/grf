/** \file optimizable.hpp 
    \brief Abstract class for optimizable objects */
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

#ifndef __OPTIMIZABLE_HPP__
#define __OPTIMIZABLE_HPP__

#include "specialtypes.hpp"

namespace bayesopt {

  
  class RBOptimizable
  {
  public:
    RBOptimizable(){};
    virtual ~RBOptimizable(){};

    virtual double evaluate(const vectord& query) = 0;
  };

  class RGBOptimizable
  {
  public:
    RGBOptimizable(){};
    virtual ~RGBOptimizable(){};

    virtual double evaluate(const vectord& query) = 0;
    virtual double evaluate(const vectord& query, 
			    vectord& grad) = 0;

  };


  class RBOptimizableWrapper
  {
  public:
    explicit RBOptimizableWrapper(RBOptimizable* rbo): rbo_(rbo){};
    virtual ~RBOptimizableWrapper(){};
    double evaluate(const vectord& query){return rbo_->evaluate(query);}
  private:
    RBOptimizable* rbo_;
  };

  class RGBOptimizableWrapper
  {
  public:
    explicit RGBOptimizableWrapper(RGBOptimizable* rgbo): rgbo_(rgbo){};
    virtual ~RGBOptimizableWrapper(){};
    double evaluate(const vectord& query, vectord& grad){return rgbo_->evaluate(query,grad);}
  private:
    RGBOptimizable* rgbo_;
  };

  
}  // namespace bayesopt

#endif
