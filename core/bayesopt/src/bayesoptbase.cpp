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
#include "bayesopt/bayesoptbase.hpp"
#include "bayesopt/parameters.hpp"

#include "log.hpp"
#include "posteriormodel.hpp"
#include "specialtypes.hpp"
#include "bopt_state.hpp"


namespace bayesopt
{
  BayesOptBase::BayesOptBase(size_t dim, Parameters parameters):
    mParameters(parameters), mDims(dim)
  {
    // Random seed
    if (mParameters.random_seed < 0) mParameters.random_seed = std::time(0); 
    mEngine.seed(mParameters.random_seed);
    
    // Setting verbose stuff (files, levels, etc.)
    int verbose = mParameters.verbose_level;
    if (verbose>=3)
      {
	FILE* log_fd = fopen( mParameters.log_filename.c_str() , "w" );
	Output2FILE::Stream() = log_fd; 
	verbose -= 3;
      }

    switch(verbose)
      {
      case 0: FILELog::ReportingLevel() = logWARNING; break;
      case 1: FILELog::ReportingLevel() = logINFO; break;
      case 2: FILELog::ReportingLevel() = logDEBUG4; break;
      default:
	FILELog::ReportingLevel() = logERROR; break;
      }
  }

  BayesOptBase::~BayesOptBase()
  { } // Default destructor


  // OPTIMIZATION INTERFACE
  void BayesOptBase::optimize(vectord &bestPoint)
  {
    assert(mDims == bestPoint.size());
    
    // Restore state from file
    if(mParameters.load_save_flag == 1 || mParameters.load_save_flag == 3)
      {
        BOptState state;
        bool load_succeed = state.loadFromFile(mParameters.load_filename, 
					       mParameters);
        if(load_succeed)
	  {
            restoreOptimization(state);
            FILE_LOG(logINFO) << "State succesfully restored from file \"" 
			      << mParameters.load_filename << "\"";
	  }
        else
	  {
	    // If load not succeed, print info message and then
	    // initialize a new optimization
            FILE_LOG(logINFO) << "File \"" << mParameters.load_filename 
			      << "\" does not exist,"
			      << " starting a new optimization";
            initializeOptimization();
	  }
      }
    else
      {
	// Initialize a new state
        initializeOptimization();
      }
    
    for (size_t ii = mCurrentIter; ii < mParameters.n_iterations; ++ii)
      {      
        stepOptimization();
      }
   
    bestPoint = getFinalResult();
  } // optimize


  void BayesOptBase::stepOptimization()
  {
    // Find what is the next point.
    vectord xNext = nextPoint(); 
    double yNext = evaluateSampleInternal(xNext);

    // If we are stuck in the same point for several iterations, try a random jump!
    if (mParameters.force_jump)
      {
        if (std::pow(mYPrev - yNext,2) < mParameters.noise)
          {
            mCounterStuck++;
            FILE_LOG(logDEBUG) << "Stuck for "<< mCounterStuck << " steps";
          }
        else
          {
            mCounterStuck = 0;
          }
        mYPrev = yNext;

        if (mCounterStuck > mParameters.force_jump)
          {
            FILE_LOG(logINFO) << "Forced random query!";
            xNext = samplePoint();
            yNext = evaluateSampleInternal(xNext);
            mCounterStuck = 0;
          }
      }

    mModel->addSample(xNext,yNext);

    // Update surrogate model
    bool retrain = ((mParameters.n_iter_relearn > 0) && 
		    ((mCurrentIter + 1) % mParameters.n_iter_relearn == 0));

    if (retrain)  // Full update
      {
        mModel->updateHyperParameters();
        mModel->fitSurrogateModel();
      }
    else          // Incremental update
      {
        mModel->updateSurrogateModel();
      } 

    plotStepData(mCurrentIter,xNext,yNext);
    mModel->updateCriteria(xNext);
    mCurrentIter++;
    
    // Save state if required
    if(mParameters.load_save_flag == 2 || mParameters.load_save_flag == 3)
      {
        BOptState state;
        saveOptimization(state);
        state.saveToFile(mParameters.save_filename);
      }
  }
  

  void BayesOptBase::initializeOptimization()
  {
    // Posterior surrogate model
    mModel.reset(PosteriorModel::create(mDims,mParameters,mEngine));
    
    // Configure iteration parameters
    if (mParameters.n_init_samples <= 0)
      {
        mParameters.n_init_samples = 
          static_cast<size_t>(ceil(0.1*mParameters.n_iterations));	
      }
    
    size_t nSamples = mParameters.n_init_samples;

    // Generate xPoints for initial sampling
    matrixd xPoints(nSamples,mDims);
    vectord yPoints(nSamples,0);

    // Save generated xPoints before its evaluation
    generateInitialPoints(xPoints);
    saveInitialSamples(xPoints);
    mModel->setSamples(xPoints);
    
    // Save on each evaluation for safety reasons
    for(size_t i=0; i<yPoints.size(); i++)
      {
        yPoints[i] = evaluateSampleInternal(row(xPoints,i));
	//We clear the vector in the first iteration
        saveResponse(yPoints[i], i==0);
      }
    
    // Put samples into model
    mModel->setSamples(yPoints);
 
    if(mParameters.verbose_level > 0)
      {
        mModel->plotDataset(logDEBUG);
      }
    
    mModel->updateHyperParameters();
    mModel->fitSurrogateModel();
    mCurrentIter = 0;

    mCounterStuck = 0;
    mYPrev = 0.0;
  }

  vectord BayesOptBase::getFinalResult()
  {
    return remapPoint(getPointAtMinimum());
  }


  // SAVE-RESTORE INTERFACE
  void BayesOptBase::saveOptimization(BOptState &state)
  {   
    // BayesOptBase members
    state.mCurrentIter = mCurrentIter;
    state.mCounterStuck = mCounterStuck;
    state.mYPrev = mYPrev;

    state.mParameters = mParameters;

    // Samples
    state.mX = mModel->getData()->mX;
    state.mY = mModel->getData()->mY;
  }

  void BayesOptBase::restoreOptimization(BOptState state)
  {
    // Restore parameters
    mParameters = state.mParameters; 
    
    // Posterior surrogate model
    mModel.reset(PosteriorModel::create(mDims, mParameters, mEngine));
    
    // Load samples, putting mX vecOfvec into a matrixd
    matrixd xPoints(state.mX.size(),state.mX[0].size());
    vectord yPoints(state.mX.size(),0);
    for(size_t i=0; i<state.mX.size(); i++)
      {
        row(xPoints, i) = state.mX[i];
        if(i < state.mY.size())
	  {
            yPoints[i] = state.mY[i];
	  }
	else
	  {
	    // Generate remaining initial samples saving in each evaluation	    
	    yPoints[i] = evaluateSampleInternal(row(xPoints,i));
	    saveResponse(yPoints[i], false);
	  }
      }
    
    // Set loaded and generated samples
    mModel->setSamples(xPoints,yPoints);
        
    if(mParameters.verbose_level > 0)
    {
        mModel->plotDataset(logDEBUG);
    }
    
    // Calculate the posterior model
    mModel->updateHyperParameters();
    mModel->fitSurrogateModel();
    
    mCurrentIter = state.mCurrentIter;
    mCounterStuck = state.mCounterStuck;
    mYPrev = state.mYPrev;
    
    // Check if optimization has already finished
    if(mCurrentIter >= mParameters.n_iterations)
      {
        FILE_LOG(logINFO) << "Optimization has already finished, delete \"" 
			  << mParameters.load_filename 
			  << "\" or give more n_iterations in parameters."; 
      }
  }

  
  // GETTERS AND SETTERS
  // Potential inline functions. Moved here to simplify API and header
  // structure.
  ProbabilityDistribution* BayesOptBase::getPrediction(const vectord& query)
  { return mModel->getPrediction(query); };
  
  const Dataset* BayesOptBase::getData()
  { return mModel->getData(); };

  Parameters* BayesOptBase::getParameters() 
  {return &mParameters;};

  double BayesOptBase::getValueAtMinimum()
  { return mModel->getValueAtMinimum(); };

  double BayesOptBase::evaluateCriteria(const vectord& query)
  {
    if (checkReachability(query)) return mModel->evaluateCriteria(query);
    else return 0.0;
  }

  size_t BayesOptBase::getCurrentIter()
  {return mCurrentIter;};

  

  // PROTECTED
  vectord BayesOptBase::getPointAtMinimum() 
  { return mModel->getPointAtMinimum(); };

  double BayesOptBase::evaluateSampleInternal( const vectord &query )
  { 
    const double yNext = evaluateSample(remapPoint(query)); 
    if (yNext == HUGE_VAL)
      {
	throw std::runtime_error("Function evaluation out of range");
      }
    return yNext;
  }; 



  
  
  void BayesOptBase::plotStepData(size_t iteration, const vectord& xNext,
				     double yNext)
  {
    if(mParameters.verbose_level >0)
      { 
	FILE_LOG(logINFO) << "Iteration: " << iteration+1 << " of " 
			  << mParameters.n_iterations << " | Total samples: " 
			  << iteration+1+mParameters.n_init_samples ;
	FILE_LOG(logINFO) << "Query: "         << remapPoint(xNext); ;
	FILE_LOG(logINFO) << "Query outcome: " << yNext ;
	FILE_LOG(logINFO) << "Best query: "    << getFinalResult(); 
	FILE_LOG(logINFO) << "Best outcome: "  << getValueAtMinimum();
      }
  } //plotStepData


  void BayesOptBase::saveInitialSamples(matrixd xPoints)
  {
    // Save state if required
    if(mParameters.load_save_flag == 2 || mParameters.load_save_flag == 3)
      {
        BOptState state;
        saveOptimization(state);
        
        // Overwrite the state with initial samples so far
        state.mX.clear();
        for(size_t i=0; i<xPoints.size1(); i++)
	  {
            state.mX.push_back(row(xPoints,i));
	  }
        state.saveToFile(mParameters.save_filename);
      }
  }


  void BayesOptBase::saveResponse(double yPoint, bool clear)
  {
    // Save state if required
    if(mParameters.load_save_flag == 2 || mParameters.load_save_flag == 3)
      {
        BOptState state;
        saveOptimization(state);
	if (clear)
	  {
	    state.mY.clear();
	  }
	utils::append(state.mY,yPoint);
        state.saveToFile(mParameters.save_filename);
      }
  }






  // PRIVATE MEMBERS
  vectord BayesOptBase::nextPoint()
  {
    //Epsilon-Greedy exploration (see Bull 2011)
    if ((mParameters.epsilon > 0.0) && (mParameters.epsilon < 1.0))
      {
	randFloat drawSample(mEngine,realUniformDist(0,1));
	double result = drawSample();
	FILE_LOG(logINFO) << "Trying random jump with prob:" << result;
	if (mParameters.epsilon > result)
	  {
	    FILE_LOG(logINFO) << "Epsilon-greedy random query!";
	    return samplePoint();
	  }
      }

    vectord Xnext(mDims);    

    // GP-Hedge and related algorithms
    if (mModel->criteriaRequiresComparison())
      {
	bool changed = true;

	mModel->setFirstCriterium();
	while (changed)
	  {
	    findOptimal(Xnext);
	    changed = mModel->setNextCriterium(Xnext);
	  }
	std::string name = mModel->getBestCriteria(Xnext);
	FILE_LOG(logINFO) << name << " was selected.";
      }
    else  // Standard "Bayesian optimization"
      {
	FILE_LOG(logDEBUG) << "------ Optimizing criteria ------";
	findOptimal(Xnext);
      }
    return Xnext;
  }



} //namespace bayesopt

