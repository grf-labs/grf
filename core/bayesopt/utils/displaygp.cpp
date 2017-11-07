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

#include "displaygp.hpp"
#include "dataset.hpp"
#include "prob_distribution.hpp"

namespace bayesopt
{
  namespace utils
  {      
    DisplayProblem1D::DisplayProblem1D(): MatPlot()
    {
      status = NOT_READY;
    }

    void DisplayProblem1D::init(BayesOptBase* bopt, size_t dim)
    {
      if (dim != 1) 
	{ 
	  throw std::invalid_argument("Display only works for 1D problems"); 
	}

      bopt_model = bopt;
      bopt->initializeOptimization();
      size_t n_points = bopt->getData()->getNSamples();
      for (size_t i = 0; i<n_points;++i)
	{
	  const double res = bopt->getData()->getSampleY(i);
	  const vectord last = bopt->getData()->getSampleX(i);
	  ly.push_back(res);
	  lx.push_back(last(0));
	}
      state_ii = 0;    
      status = STOP;
    };

    void DisplayProblem1D::setSTEP()
    {
      if (status != NOT_READY)
	{
	  status = STEP;
	}
    };

    void DisplayProblem1D::toogleRUN()
    {
      if (status != NOT_READY)
	{
	  if(status != RUN)
	    {
	      status = RUN;
	    }
	  else
	    {
	      status = STOP;
	    }
	}
    }

    void DisplayProblem1D::DISPLAY()
    {
      if (status != NOT_READY)
	{
	  size_t nruns = bopt_model->getParameters()->n_iterations;
	  if ((status != STOP) && (state_ii < nruns))
	    {
	      // We are moving. Next iteration
	      ++state_ii;
	      bopt_model->stepOptimization(); 
	      const double res = bopt_model->getData()->getLastSampleY();
	      const vectord last = bopt_model->getData()->getLastSampleX();
	      ly.push_back(res);
	      lx.push_back(last(0));
	      
	      if (status == STEP) { status = STOP; }
	    }
	  
	  // We compute the prediction, true value and criteria at 1000 points
	  int n=1000;
	  std::vector<double> x,y,z,su,sl,c;
	  x = linspace(0,1,n);
	  y = x; z = x; su = x; sl = x; c = x;
	  
	  // Query functions at the 1000 points
	  vectord q(1);
	  for(size_t i=0; i<n; ++i)
	    {
	      q(0) = x[i];                                                 // Query
	      ProbabilityDistribution* pd = bopt_model->getPrediction(q);
	      y[i] = pd->getMean();                                //Expected value
	      su[i] = y[i] + 2*pd->getStd();                       //Upper bound (95 %)
	      sl[i] = y[i] - 2*pd->getStd();                       //Lower bound (95 %)
	      c[i] = -bopt_model->evaluateCriteria(q);             //Criteria value
	      z[i] = bopt_model->evaluateSample(q);                //Target function true value
	    }
 
	  //GP subplot
	  subplot(2,1,1);
	  title("Press r to run and stop, s to run a step and q to quit.");
	  plot(x,y); set(3);                            // Expected value in default color (blue)
	  plot(lx,ly);set("k");set("o");set(4);         // Data points as black star
	  plot(x,su);set("g"); set(2);                  // Uncertainty as green lines
	  plot(x,sl);set("g"); set(2);
	  plot(x,z);set("r"); set(3);                   // True function as red line
	  
	  //Criterion subplot
	  subplot(2,1,2);
	  plot(x,c); set(3);
	}
    };

    DisplayProblem2D::DisplayProblem2D(): 
      MatPlot(), cx(1), cy(1), c_points(100), cX(c_points),
      cY(c_points), cZ(c_points,std::vector<double>(c_points))
    {
      status = NOT_READY;
    }

    void DisplayProblem2D::setSolution(vectord sol)
    {
      solx.push_back(sol(0));
      soly.push_back(sol(1));
    }

    void DisplayProblem2D::prepareContourPlot()
    {
      cX=linspace(0,1,c_points);
      cY=linspace(0,1,c_points);
      
      for(int i=0;i<c_points;++i)
	{
	  for(int j=0;j<c_points;++j)
	    {
	      vectord q(2);
	      q(0) = cX[j]; q(1) = cY[i];
	      cZ[i][j]= bopt_model->evaluateSample(q);
	    }
	}
    }

    void DisplayProblem2D::init(BayesOptBase* bopt, size_t dim)
    {
      if (dim != 2) 
	{ 
	  throw std::invalid_argument("This display only works "
				      "for 2D problems"); 
	}
      
      bopt_model = bopt;
      prepareContourPlot();
      
      bopt->initializeOptimization();
      size_t n_points = bopt->getData()->getNSamples();
      for (size_t i = 0; i<n_points;++i)
	{
	  const vectord last = bopt->getData()->getSampleX(i);
	  lx.push_back(last(0));
	  ly.push_back(last(1));
	}
      state_ii = 0;    
      status = STOP;
    };

    void DisplayProblem2D::setSTEP()
    {
      if (status != NOT_READY)
	{
	  status = STEP;
	}
    };

    void DisplayProblem2D::toogleRUN()
    {
      if (status != NOT_READY)
	{
	  if(status != RUN)
	    {
	      status = RUN;
	    }
	  else
	    {
	      status = STOP;
	    }
	}
    }
    
    void DisplayProblem2D::DISPLAY()
    {
      if (status != NOT_READY)
	{
	  size_t nruns = bopt_model->getParameters()->n_iterations;
	  title("Press r to run and stop, s to run a step and q to quit.");
	  contour(cX,cY,cZ,50);                         // Contour plot (50 lines)
	  plot(cx,cy);set("g");set("o");set(4);         // Data points as black star
	  plot(solx,soly);set("r"); set("o");set(4);    // Solutions as red points
	  
	  if ((status != STOP) && (state_ii < nruns))
	    {
	      // We are moving. Next iteration
	      ++state_ii;
	      bopt_model->stepOptimization(); 
	      const vectord last = bopt_model->getData()->getLastSampleX();
	      //GP subplot
	      cx[0] = last(0);
	      cy[0] = last(1);
	      
	      if (!lx.empty())
		{	
		  plot(lx,ly);set("k");set("o");set(4);         // Data points as black star
		}
	      
	      lx.push_back(last(0));
	      ly.push_back(last(1));
	      
	      if (status == STEP) { status = STOP; }
	    }	    
	  else
	    {
	      plot(lx,ly);set("k");set("o");set(4);         // Data points as black star
	    }
	  
	}
    };

  } //namespace utils

} //namespace bayesopt
