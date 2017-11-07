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

int main(int nargs, char *args[])
{
  bayesopt::Parameters par;
  par.verbose_level = 0;
  par.noise = 1e-10;
  par.force_jump = 30;

  std::ofstream log;
  std::clock_t start_t;


  /* Branin */
  log.open("branin.log");
  par.n_init_samples = 5;
  par.n_iterations = 195;

  for (size_t ii = 0; ii < 10; ++ii)
    {
      par.random_seed = ii;
      BraninNormalized branin(par);
      vectord result(2);

      start_t = clock();
      branin.initializeOptimization();
      
      for (size_t jj = 0; jj < par.n_iterations; ++jj)
  	{      
  	  branin.stepOptimization();
  	  if (jj == 50)
  	    {
  	      result = branin.getFinalResult();	      
  	      log << branin.evaluateSample(result) << ", ";
  	    }
  	}
      result = branin.getFinalResult();	      
      log << branin.evaluateSample(result) << ", ";
      
      log << static_cast<double>(clock() - start_t) / static_cast<double>(CLOCKS_PER_SEC)
  	  << std::endl;
      }

  log.close();


  /* Camel */
  log.open("camel.log");
  par.n_init_samples = 5;
  par.n_iterations = 95;

  for (size_t ii = 0; ii < 10; ++ii)
    {
      par.random_seed = ii;
      ExampleCamelback camel(par);
      vectord result(2);

      vectord lb(2); lb(0) = -2; lb(1) = -1;
      vectord ub(2); ub(0) =  2; ub(1) = 1;

      camel.setBoundingBox(lb,ub);

      start_t = clock();
      camel.initializeOptimization();
      
      for (size_t jj = 0; jj < par.n_iterations; ++jj)
  	{      
  	  camel.stepOptimization();
  	  if (jj == 50)
  	    {
  	      result = camel.getFinalResult();	      
  	      log << camel.evaluateSample(result) << ", ";
  	    }
  	}
      result = camel.getFinalResult();	      
      log << camel.evaluateSample(result) << ", ";
      
      log << static_cast<double>(clock() - start_t) / static_cast<double>(CLOCKS_PER_SEC)
  	  << std::endl;
      }

  log.close();


  /* Hart */
  log.open("hart.log");
  par.n_init_samples = 10;
  par.n_iterations = 190;

  for (size_t ii = 0; ii < 10; ++ii)
    {
      par.random_seed = ii;
      ExampleHartmann6 hart(par);
      vectord result(6);

      start_t = clock();
      hart.initializeOptimization();
      
      for (size_t jj = 0; jj < par.n_iterations; ++jj)
  	{      
  	  hart.stepOptimization();
  	  if (jj == 50)
  	    {
  	      result = hart.getFinalResult();	      
  	      log << hart.evaluateSample(result) << ", ";
  	    }
  	}
      result = hart.getFinalResult();	      
      log << hart.evaluateSample(result) << ", ";
      
      log << static_cast<double>(clock() - start_t) / static_cast<double>(CLOCKS_PER_SEC)
  	  << std::endl;
      }

  log.close();


  /***********************************************************************/
  par.n_init_samples = 2;
  par.n_iter_relearn = 1;
  
  par.l_type = L_MCMC;
  par.sc_type = SC_MAP;


  /* Branin */
  log.open("branin_mcmc.log");
  par.n_iterations = 198;

  for (size_t ii = 0; ii < 10; ++ii)
    {
      par.random_seed = ii;
      BraninNormalized branin(par);
      vectord result(2);

      start_t = clock();
      branin.initializeOptimization();
      
      for (size_t jj = 0; jj < par.n_iterations; ++jj)
	{      
	  branin.stepOptimization();
	  if (jj == 50)
	    {
	      result = branin.getFinalResult();	      
	      log << branin.evaluateSample(result) << ", ";
	    }
	}
      result = branin.getFinalResult();	      
      log << branin.evaluateSample(result) << ", ";
      
      log << static_cast<double>(clock() - start_t) / static_cast<double>(CLOCKS_PER_SEC)
	  << std::endl;
      }

  log.close();


  /* Camel */
  log.open("camel_mcmc.log");
  par.n_iterations = 98;

  for (size_t ii = 0; ii < 10; ++ii)
    {
      par.random_seed = ii;
      ExampleCamelback camel(par);
      vectord result(2);

      vectord lb(2); lb(0) = -2; lb(1) = -1;
      vectord ub(2); ub(0) =  2; ub(1) = 1;

      camel.setBoundingBox(lb,ub);

      start_t = clock();
      camel.initializeOptimization();
      
      for (size_t jj = 0; jj < par.n_iterations; ++jj)
	{      
	  camel.stepOptimization();
	  if (jj == 50)
	    {
	      result = camel.getFinalResult();	      
	      log << camel.evaluateSample(result) << ", ";
	    }
	}
      result = camel.getFinalResult();	      
      log << camel.evaluateSample(result) << ", ";
      
      log << static_cast<double>(clock() - start_t) / static_cast<double>(CLOCKS_PER_SEC)
	  << std::endl;
      }

  log.close();


  /* Hart */
  log.open("hart_mcmc.log");
  par.n_iterations = 198;

  for (size_t ii = 0; ii < 10; ++ii)
    {
      par.random_seed = ii;
      ExampleHartmann6 hart(par);
      vectord result(6);

      start_t = clock();
      hart.initializeOptimization();
      
      for (size_t jj = 0; jj < par.n_iterations; ++jj)
	{      
	  hart.stepOptimization();
	  if (jj == 50)
	    {
	      result = hart.getFinalResult();	      
	      log << hart.evaluateSample(result) << ", ";
	    }
	}
      result = hart.getFinalResult();	      
      log << hart.evaluateSample(result) << ", ";
      
      log << static_cast<double>(clock() - start_t) / static_cast<double>(CLOCKS_PER_SEC)
	  << std::endl;
      }

  log.close();


  return 0;
}
