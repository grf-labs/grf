/*-------------------------------------------------------------------------------
 This file is part of Ranger.

 Ranger is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 Ranger is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with Ranger. If not, see <http://www.gnu.org/licenses/>.

 Written by:

 Marvin N. Wright
 Institut für Medizinische Biometrie und Statistik
 Universität zu Lübeck
 Ratzeburger Allee 160
 23562 Lübeck

 http://www.imbs-luebeck.de
 wright@imbs.uni-luebeck.de
 #-------------------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <utility/DataDouble.h>
#include <utility/DataChar.h>
#include <utility/DataFloat.h>

#include "ArgumentHandler.h"
#include "relabeling/QuantileRelabelingStrategy.h"
#include "relabeling/InstrumentalRelabelingStrategy.h"
#include "prediction/QuantilePredictionStrategy.h"
#include "prediction/InstrumentalPredictionStrategy.h"
#include "splitting/RegressionSplittingRule.h"
#include "splitting/ProbabilitySplittingRule.h"
#include "ForestModel.h"

Data* initializeData(ArgumentHandler& arg_handler) {
  Data* data = NULL;

  // Initialize data with memmode
  switch (arg_handler.memmode) {
    case MEM_DOUBLE:
      data = new DataDouble();
      break;
    case MEM_FLOAT:
      data = new DataFloat();
      break;
    case MEM_CHAR:
      data = new DataChar();
      break;
  }

  bool rounding_error = data->loadFromFile(arg_handler.file);
  if (rounding_error) {
    std::cerr << "Warning: Rounding or Integer overflow occurred."
        "Use FLOAT or DOUBLE precision to avoid this." << std::endl;
  }
  data->sort();

  return data;
}

//int main(int argc, char **argv) {
//
//  ArgumentHandler arg_handler(argc, argv);
//  ForestModel* forest_model = 0;
//  try {
//
//    // Handle command line arguments
//    if (arg_handler.processArguments() != 0) {
//      return 0;
//    }
//    arg_handler.checkArguments();
//
//    Data* data = initializeData(arg_handler);
//    size_t dependent_varID = data->getVariableID(arg_handler.depvarname);
//
//    switch (arg_handler.treetype) {
//    case TREE_QUANTILE: {
//      std::vector<double> *default_quantiles = new std::vector<double>({0.15, 0.5, 0.85});
//      std::vector<double> *quantiles = !arg_handler.quantiles->empty()
//                                       ? arg_handler.quantiles
//                                       : default_quantiles;
//
//      std::unordered_map<std::string, size_t> observables = {{"outcome", dependent_varID}};
//      RelabelingStrategy* relabeling_strategy = new QuantileRelabelingStrategy(quantiles, dependent_varID);
//      SplittingRule* splitting_rule = new ProbabilitySplittingRule(data, quantiles->size());
//      PredictionStrategy* prediction_strategy = new QuantilePredictionStrategy(quantiles);
//
//      forest_model = new ForestModel(observables,
//                                     relabeling_strategy,
//                                     splitting_rule,
//                                     prediction_strategy);
//      break;
//    }
//    case TREE_INSTRUMENTAL: {
//      size_t treatment_varID = 0;
//      size_t instrument_varID = 0;
//      if (arg_handler.predict.empty()) {
//        treatment_varID = data->getVariableID(arg_handler.statusvarname);
//        instrument_varID = data->getVariableID(arg_handler.instrumentvarname);
//      }
//
//      std::unordered_map<std::string, size_t> observables = {
//          {"outcome", dependent_varID},
//          {"treatment", treatment_varID},
//          {"instrument", instrument_varID}};
//
//      RelabelingStrategy *relabeling_strategy = new InstrumentalRelabelingStrategy(observables);
//      SplittingRule *splitting_rule = new RegressionSplittingRule(data);
//      PredictionStrategy* prediction_strategy = new InstrumentalPredictionStrategy();
//
//      forest_model = new ForestModel(observables,
//                                     relabeling_strategy,
//                                     splitting_rule,
//                                     prediction_strategy);
//      break;
//    }
//      default:
//        throw std::runtime_error("Unrecognized tree type.");
//    }
//
//
//    // Verbose output to logfile if non-verbose mode
//    std::ostream* verbose_out;
//    if (arg_handler.verbose) {
//      verbose_out = &std::cout;
//    } else {
//      std::ofstream* logfile = new std::ofstream();
//      logfile->open("ranger.log");
//      if (!logfile->good()) {
//        throw std::runtime_error("Could not write to logfile.");
//      }
//      verbose_out = logfile;
//    }
//
//    // Call Ranger
//    *verbose_out << "Starting Ranger." << std::endl;
//    forest_model->initCpp(arg_handler.memmode, arg_handler.mtry,
//        arg_handler.ntree, verbose_out, arg_handler.seed, arg_handler.nthreads,
//        arg_handler.predict, arg_handler.targetpartitionsize, arg_handler.splitweights,
//        arg_handler.alwayssplitvars, arg_handler.replace,
//        arg_handler.savemem, arg_handler.caseweights, arg_handler.fraction);
//
//    Forest* forest = forest_model->train(data);
//    std::vector<std::vector<double>> predictions = forest_model->predict(forest, data);
//    forest_model->writeOutput(data, predictions);
//    *verbose_out << "Finished Ranger." << std::endl;
//
//    delete forest_model;
//    delete forest;
//  } catch (std::exception& e) {
//    std::cerr << "Error: " << e.what() << " Ranger will EXIT now." << std::endl;
//    delete forest_model;
//    return -1;
//  }
//
//  return 0;
//}
