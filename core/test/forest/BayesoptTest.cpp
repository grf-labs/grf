/*-------------------------------------------------------------------------------
  This file is part of generalized random forest (grf).

  grf is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  grf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with grf. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#include "commons/utility.h"
#include "forest/ForestPredictor.h"

#include "forest/BayesoptSrc.h"

#include <ctime>
//#include "bayesopt.hpp"
#include <bayesopt/bayesopt.hpp>              // For the C++ API
#include <bayesopt/parameters.hpp>              // For the C++ API
#include <boost/numeric/ublas/assignment.hpp> // <<= op assigment

#include "catch.hpp"

TEST_CASE("testing auto parameter tuning with Bayesopt", "[tuning, bayesopt]") {

     //*/
    int n = 3;                   // Number of dimensions
    clock_t start, end;
    double diff,diff2;

    // Common configuration
    // See parameters.h for the available options.
    // Some parameters did not need to be changed for default, but we have done it for
    // illustrative purpose.
    bayesopt::Parameters par = initialize_parameters_to_default();

    //*
    // budget parameters
    par.n_iterations = 10;    // Number of iterations
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
    std::cout << "Running C++ interface" << std::endl;
    //*
    //Data* data = load_data("/Users/kuangkun/Documents/Stanford/Generalized Random Forests/Code/grftuning/core/test/forest/dataset/AirFoil.txt");
    Data* data = load_data("../test/forest/dataset/AirFoil.txt");
    uint outcome_index = data->get_num_cols()-1;
    BayesoptSrc opt(n,par);
    opt.setData(data, outcome_index);
    vectord result(n);

    // set the boundary for tuning parameters
    boost::numeric::ublas::vector<double> lowerBound(n);
    boost::numeric::ublas::vector<double> upperBound(n);
    lowerBound[0] = 0.01;lowerBound[1] = 1;lowerBound[2] = 0.01;
    upperBound[0] = 1.0;upperBound[1] = 100;upperBound[2] = 0.5;
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
}
