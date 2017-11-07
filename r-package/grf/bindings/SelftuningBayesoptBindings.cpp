#include <map>
#include <Rcpp.h>
#include <sstream>
#include <vector>

#include "commons/utility.h"
#include "RcppUtilities.h"
#include "forest/BayesoptSrc.h"

#include <bayesopt/bayesopt.hpp>              // For the C++ API
#include <bayesopt/parameters.hpp>              // For the C++ API
#include <boost/numeric/ublas/assignment.hpp> // <<= op assigment



// [[Rcpp::export]]
Rcpp::List selftuning_train(Rcpp::NumericMatrix input_data,
                            size_t outcome_index,
                            std::vector <std::string> variable_names,
                            double mtry_l, double mtry_r, 
                            double min_node_size_l, double min_node_size_r, 
                            double sample_faction_l, double sample_faction_r) {
  Data* data = RcppUtilities::convert_data(input_data, variable_names);
  //Data* data = load_data("/Users/kuangkun/Documents/Stanford/Generalized Random Forests/Code/grftuning_bak/core/test/forest/dataset/AirFoil.txt");
  
  outcome_index = outcome_index - 1; // from R to C++, the index need to substract 1.

  //*/
  int n = 3;                   // Number of dimensions
  clock_t start, end;
  double diff;

  // Common configuration
  // See parameters.h for the available options.
  // Some parameters did not need to be changed for default, but we have done it for
  // illustrative purpose.
  bayesopt::Parameters par = initialize_parameters_to_default();

  // budget parameters
  par.n_iterations = 100;    // Number of iterations
  par.random_seed = 10;
  par.n_init_samples = 10;
  par.n_iter_relearn = 5;
  par.init_method = 2; //"Sobol sequences"


  //Exploration/Exploitation parameters
  par.crit_name = "cEI";
  par.epsilon = 0.1;

  //Surrogate parameters
  par.surr_name = "sGaussianProcessML";
  par.kernel.name = "kMaternARD5";

  //hyperparameter learning
  par.sc_type = SC_MAP;
  par.l_type = L_MCMC;
  //*/

  /*******************************************/

  //*
  BayesoptSrc::BayesoptSrc opt(n,par);
  opt.setData(data,outcome_index);
  vectord result(n);

  // set the boundary for tuning parameters
  boost::numeric::ublas::vector<double> lowerBound(n);
  boost::numeric::ublas::vector<double> upperBound(n);
  lowerBound[0] = mtry_l;lowerBound[1] = min_node_size_l;lowerBound[2] = sample_faction_l;
  upperBound[0] = mtry_r;upperBound[1] = min_node_size_r;upperBound[2] = sample_faction_r;
  std::cout << "****" << lowerBound << upperBound << std::endl;
  opt.setBoundingBox(lowerBound, upperBound);

  // Run C++ interface
  start = clock();
  opt.optimize(result);
  end = clock();
  diff = (double)(end-start) / (double)CLOCKS_PER_SEC;


  // Results
  std::cout << "Best Query: " << result << std::endl;
  std::cout << "Best Outcome: " << opt.evaluateSample(result) << std::endl;
  //std::cout << "Best Outcome with 5000 trees: " << opt.evaluateSample_with_5000_trees(result) << std::endl;
  std::cout << "Elapsed time: " << diff << " seconds" << std::endl;
  //*/

  Rcpp::List result_return;
  result_return.push_back(result, "bestParameters");
  delete data;
  return result_return;
}
