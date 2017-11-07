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

#include "log.hpp"

#include "gaussian_process.hpp"
#include "gaussian_process_ml.hpp"
#include "gaussian_process_normal.hpp"
#include "student_t_process_jef.hpp"
#include "student_t_process_nig.hpp"

#include "nonparametricprocess.hpp"


namespace bayesopt
{

  NonParametricProcess::NonParametricProcess(size_t dim, Parameters parameters, 
					     const Dataset& data, 
					     MeanModel& mean,
					     randEngine& eng):
    mData(data), dim_(dim), mMean(mean), mSigma(parameters.sigma_s)
  {}

  NonParametricProcess::~NonParametricProcess(){}


  NonParametricProcess* NonParametricProcess::create(size_t dim, 
						     Parameters parameters, 
						     const Dataset& data, 
						     MeanModel& mean,
						     randEngine& eng)
  {
    NonParametricProcess* s_ptr;

    std::string name = parameters.surr_name;

    if (!name.compare("sGaussianProcess"))
      s_ptr = new GaussianProcess(dim,parameters,data,mean,eng);
    else  if(!name.compare("sGaussianProcessML"))
      s_ptr = new GaussianProcessML(dim,parameters,data,mean,eng);
    else  if(!name.compare("sGaussianProcessNormal"))
      s_ptr = new GaussianProcessNormal(dim,parameters,data,mean,eng);
    else if (!name.compare("sStudentTProcessJef"))
      s_ptr = new StudentTProcessJeffreys(dim,parameters,data,mean,eng); 
    else if (!name.compare("sStudentTProcessNIG"))
      s_ptr = new StudentTProcessNIG(dim,parameters,data,mean,eng); 
    else
      {
	throw std::invalid_argument("Surrogate function not supported");
      }
    return s_ptr;
  };

} //namespace bayesopt
