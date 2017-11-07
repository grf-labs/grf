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
#include "param_loader.hpp"

int main(int nargs, char *args[])
{
  bayesopt::Parameters parameters;
  if(nargs > 1){
    if(!bayesopt::utils::ParamLoader::load(args[1], parameters)){
        std::cout << "ERROR: provided file \"" << args[1] << "\" does not exist" << std::endl;
        return -1;
    }
  }
  else{
    parameters.n_init_samples = 10;
    parameters.n_iterations = 300;

    parameters.crit_name = "cHedge(cEI,cLCB,cExpReturn,cOptimisticSampling)";
    //bayesopt::utils::ParamLoader::save("bo_oned.txt", parameters);
  }
  
  ExampleOneD opt(parameters);
  vectord result(1);
  opt.optimize(result);
  
  std::cout << "Result:" << result << std::endl;
  opt.printOptimal();

  return 0;
}
