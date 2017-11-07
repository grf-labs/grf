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

#include <ctime>
#include "bayesopt/bayesopt.h"                 // For the C API
#include "bayesopt/bayesopt.hpp"               // For the C++ API
#include <boost/numeric/ublas/assignment.hpp> // <<= op assigment


/* Function to be used for C-API testing */
double testFunction(unsigned int n, const double *x,
        double *gradient, /* NULL if not needed */
        void *func_data)
{
  double f = 10.;
  for (unsigned int i = 0; i < n; ++i)
    {
      f += (x[i] - .53) * (x[i] - .53);
    }
  return f;
}

/* Class to be used for C++-API testing */
class ExampleQuadratic: public bayesopt::ContinuousModel
{
 public:

  ExampleQuadratic(size_t dim,bayesopt::Parameters param):
    ContinuousModel(dim,param) {}

  double evaluateSample( const vectord &Xi ) 
  {
    double x[100];
    for (size_t i = 0; i < Xi.size(); ++i)
      {
  x[i] = Xi(i); 
      }
    return testFunction(Xi.size(),x,NULL,NULL);
  };


  bool checkReachability( const vectord &query )
  { return true; };
 
};


int main(int nargs, char *args[])
{    
  int n = 10;                   // Number of dimensions
  clock_t start, end;
  double diff,diff2;

  // Common configuration
  // See parameters.h for the available options.
  // Some parameters did not need to be changed for default, but we have done it for
  // illustrative purpose.
  bayesopt::Parameters par = initialize_parameters_to_default();

  par.kernel.name = "kSum(kSEISO,kConst)";
  par.kernel.hp_mean <<= 1.0, 1.0;
  par.kernel.hp_std <<= 1.0, 1.0;

  par.mean.name = "mConst";
  par.mean.coef_mean <<= 1.0;
  par.mean.coef_std <<= 1.0;
  

  par.surr_name = "sStudentTProcessJef";
  par.noise = 1e-10;

  par.sc_type = SC_MAP;
  par.l_type = L_EMPIRICAL;

  par.n_iterations = 100;    // Number of iterations
  par.random_seed = 0;
  par.n_init_samples = 15;
  par.n_iter_relearn = 0;

  /*******************************************/
  std::cout << "Running C++ interface" << std::endl;

  ExampleQuadratic opt(n,par);
  vectord result(n);

  // Run C++ interface
  start = clock();
  opt.optimize(result);
  end = clock();
  diff = (double)(end-start) / (double)CLOCKS_PER_SEC;

  /*******************************************/
  std::cout << "Running C inferface" << std::endl;
  
  // Prepare C interface
  double low[128], up[128], xmin[128], fmin[128];

  // Lower and upper bounds
  for (int i = 0; i < n; ++i) 
    {
      low[i] = 0.;    
      up[i] = 1.;
    }

  // Run C interface
  start = clock();
  bayes_optimization(n,&testFunction,NULL,low,up,xmin,fmin,par.generate_bopt_params());
  end = clock();
  diff2 = (double)(end-start) / (double)CLOCKS_PER_SEC;
  /*******************************************/


  // Results
  std::cout << "Final result C++: " << result << std::endl;
  std::cout << "Elapsed time in C++: " << diff << " seconds" << std::endl;

  std::cout << "Final result C: [" << n <<"](" << xmin[0];
  for (int i = 1; i < n; ++i )
    {
      std::cout << "," << xmin[i];      
    }
  std::cout << ")" << std::endl;
  std::cout << "Elapsed time in C: " << diff2 << " seconds" << std::endl;

}