/** \file gaussian_process_hierarchical.hpp 
    \brief Hierarchical model for Gaussian process. Abstract module */
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


#ifndef __HIERARCHICAL_GAUSSIAN_PROCESS_HPP__
#define __HIERARCHICAL_GAUSSIAN_PROCESS_HPP__

#include "conditionalbayesprocess.hpp"


namespace bayesopt
{
  
  /** \addtogroup NonParametricProcesses */
  /**@{*/

  /**
   * \brief Virtual class for hierarchical Gaussian processes.
   */
  class HierarchicalGaussianProcess: public ConditionalBayesProcess
  {
  public:
    HierarchicalGaussianProcess(size_t dim, Parameters params, const Dataset& data, 
				MeanModel& mean,randEngine& eng);
    virtual ~HierarchicalGaussianProcess() {};

  protected:
    /** 
     * \brief Computes the negative log likelihood of the data for all
     * the parameters.
     * @return value negative log likelihood
     */
    double negativeTotalLogLikelihood();

  };

  /**@}*/

} //namespace bayesopt

#endif
