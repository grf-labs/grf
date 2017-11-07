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

#include "testfunctions.hpp"
#include <ctime>
#include <fstream>
#include "param_loader.hpp"

int main(int nargs, char *args[])
{
  bayesopt::Parameters par;
  if(nargs > 1){
    if(!bayesopt::utils::ParamLoader::load(args[1], par)){
        std::cout << "ERROR: provided file \"" << args[1] << "\" does not exist" << std::endl;
        return -1;
    }
  }
  else{
    par = initialize_parameters_to_default();
    par.n_iterations = 190;
    //  par.n_iter_relearn = 0;
    par.random_seed = 0;
    par.verbose_level = 1;
    par.noise = 1e-10;
    //bayesopt::utils::ParamLoader::save("bo_branin_timed.txt", par);  
  }

  BraninNormalized branin(par);

  std::ofstream timelog;
  timelog.open("time_branin.log");
  std::clock_t curr_t;
  std::clock_t prev_t = clock();

  branin.initializeOptimization();
      
  for (size_t ii = 0; ii < par.n_iterations; ++ii)
    {      
      branin.stepOptimization();

      curr_t = clock();
      timelog << ii << ","
	      << static_cast<double>(curr_t - prev_t) / CLOCKS_PER_SEC 
	      << std::endl;
      prev_t = curr_t;
      }

  timelog.close();

  vectord result = branin.getFinalResult();
  std::cout << "Result: " << result << "->" 
	    << branin.evaluateSample(result) << std::endl;
  branin.printOptimal();

  return 0;
}
