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

#include <boost/numeric/ublas/matrix_proxy.hpp>  // For row()

#include "bayesopt/bayesopt.h"               // For the C API
#include "bayesopt/bayesopt.hpp"             // For the C++ API
#include "lhs.hpp"
#include "specialtypes.hpp"


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
class ExampleDisc: public bayesopt::DiscreteModel
{
 public:

  ExampleDisc(const vecOfvec & validSet, bayesopt::Parameters param):
    DiscreteModel(validSet,param) {}

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
  // Common configuration
  // See parameters.h for the available options
  // This is just an example of what can be done.
  bayesopt::Parameters par;

  par.surr_name = "sStudentTProcessJef";
  par.noise = 1e-10;
  par.n_iterations = 70;       // Number of iterations
  par.n_init_samples = 20;
  /*******************************************/
  
  const size_t nPoints = 1000;          // Number of discrete points
  const size_t n = 6;                   // Number of dimensions

  // Let's generate a discrete grid by latin hypercube sampling in the n-dimensional
  // space.
  randEngine mtRandom;
  matrixd xPoints(nPoints,n);

  bayesopt::utils::lhs(xPoints,mtRandom);


  // Set of discrete points using C arrays.  
  // Thanks to the quirks of Visual Studio, the following expression is invalid, so we
  // have to replace it by a literal.  WARNING: Make sure to update the size of the array
  // if the number of points or dimensions change.
#ifdef _MSC_VER
  double xPointsArray[6000];
#else
  const size_t nPinArr = n*nPoints;
  double xPointsArray[nPinArr];
#endif
  
  // Set of discrete points using ublas vectors. For the C++ interface.
  vecOfvec xP;   


  // Lets copy the generated points to the specific structure for each interface.
  for(size_t i = 0; i<nPoints; ++i)
    {
      vectord point = row(xPoints,i);  
      xP.push_back(point);
      for(size_t j=0; j<n; ++j)
	{
	  xPointsArray[i*n+j] = point(j);	  
	}
    }
    

  /*******************************************/
  // Run C++ interface
  std::cout << "Running C++ interface" << std::endl;
  ExampleDisc opt(xP,par);
  vectord result(n);
  opt.optimize(result);

  
  // Run C interface
  std::cout << "Running C interface" << std::endl;
  double x[128], fmin[128];
  bayes_optimization_disc(n, &testFunction, NULL, xPointsArray, nPoints,
			  x, fmin, par.generate_bopt_params());


  // Find the optimal value by brute force
  size_t min = 0;
  double minvalue = opt.evaluateSample(row(xPoints,min));
  for(size_t i = 1; i<nPoints; ++i)
    {
      vectord point = row(xPoints,i);  
      if (opt.evaluateSample(point) < minvalue)
	{
	  min = i;
	  minvalue = opt.evaluateSample(row(xPoints,min));
	  std::cout << i << "," << minvalue << std::endl;
	}
    }

  std::cout << "Final result C++: " << result << " | Value:" << opt.evaluateSample(result) << std::endl;
  std::cout << "Final result C: ["<< n << "](";
  for (int i = 0; i < n; i++ )
    std::cout << x[i] << ", ";
  std::cout << ")" << " | Value:" << fmin[0] << std::endl;
  std::cout << "Optimal: " << row(xPoints,min) << " | Value:" << minvalue << std::endl;
}
