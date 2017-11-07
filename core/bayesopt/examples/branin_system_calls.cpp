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

#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <boost/numeric/ublas/assignment.hpp>
#include "bayesopt/bayesopt.hpp"
#include "param_loader.hpp"
#include "fileparser.hpp"

#ifndef M_PI
/* It shouldn't be necessary, but Windows is completely nuts about
   math constants and I HATE when include order matters. */
    #define M_PI       3.14159265358979323846
#endif

class SystemCallsBranin: public bayesopt::ContinuousModel
{
public:
  SystemCallsBranin(bayesopt::Parameters par):
    ContinuousModel(2,par) {}

  double evaluateSample( const vectord& xin)
  {
     if (xin.size() != 2)
      {
	std::cout << "WARNING: This only works for 2D inputs." << std::endl
		  << "WARNING: Using only first two components." << std::endl;
      }

    float y = -1;
    double x1 = xin(0);
    double x2 = xin(1);
    
    bayesopt::utils::FileParser fp("results.txt");
    std::string call = "python ../examples/standalone_calls/eval_branin.py " + 
      fp.to_string(x1) + " " + fp.to_string(x2);
    int ret = system(call.c_str());
    
    fp.openInput();
    fp.read("y",y);
    fp.close();
    
    return y;
  }

  bool checkReachability(const vectord &query)
  {return true;};

  inline double sqr( double x ){ return x*x; };

  void printOptimal()
  {
    vectord sv(2);  
    sv(0) = 0.1238938; sv(1) = 0.818333;
    std::cout << "Solutions: " << sv << "->" 
	      << evaluateSample(sv) << std::endl;
    sv(0) = 0.5427728; sv(1) = 0.151667;
    std::cout << "Solutions: " << sv << "->" 
	      << evaluateSample(sv) << std::endl;
    sv(0) = 0.961652; sv(1) = 0.1650;
    std::cout << "Solutions: " << sv << "->" 
	      << evaluateSample(sv) << std::endl;
  }

};

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
    par.random_seed = 0;
    par.verbose_level = 1;
    par.noise = 1e-10;
    //bayesopt::utils::ParamLoader::save("system_opt.txt", par);
  }

    
  SystemCallsBranin branin(par);
  vectord result(2);

  branin.optimize(result);
  std::cout << "Result: " << result << "->" 
	    << branin.evaluateSample(result) << std::endl;
  branin.printOptimal();
  
  // Remove results.txt file
  std::string filename("results.txt");
  if( remove( filename.c_str() ) == 0 ){
    std::cout << "File \"" << filename << "\" successfully removed" << std::endl;
  }
  else{
    std::cout << "Error: cannot remove \"" << filename << "\" file" << std::endl; 
  }
  
  return 0;
}


