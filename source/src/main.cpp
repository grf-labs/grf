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
#include <Tree/QuantileRelabelingStrategy.h>
#include <Tree/ProbabilitySplittingRule.h>
#include <Tree/InstrumentalRelabelingStrategy.h>
#include <Tree/RegressionSplittingRule.h>
#include <utility/DataDouble.h>
#include <utility/DataChar.h>
#include <utility/DataFloat.h>

#include "ArgumentHandler.h"
#include "ForestQuantile.h"
#include "ForestInstrumental.h"

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

int main(int argc, char **argv) {

  ArgumentHandler arg_handler(argc, argv);
  Forest* forest = 0;
  try {

    // Handle command line arguments
    if (arg_handler.processArguments() != 0) {
      return 0;
    }
    arg_handler.checkArguments();

    Data* data = initializeData(arg_handler);
    size_t dependent_varID = data->getVariableID(arg_handler.depvarname);

    switch (arg_handler.treetype) {
    case TREE_QUANTILE: {
      std::vector<double> *default_quantiles = new std::vector<double>({0.15, 0.5, 0.85});
      std::vector<double> *quantiles = !arg_handler.quantiles->empty()
                                       ? arg_handler.quantiles
                                       : default_quantiles;

      RelabelingStrategy* relabeling_strategy = new QuantileRelabelingStrategy(quantiles, dependent_varID);
      SplittingRule* splitting_rule = new ProbabilitySplittingRule(data, quantiles->size());
      forest = new ForestQuantile(quantiles, relabeling_strategy, splitting_rule);
      break;
    }
    case TREE_INSTRUMENTAL: {
      size_t treatment_varID = 0;
      size_t instrument_varID = 0;
      if (arg_handler.predict.empty()) {
        treatment_varID = data->getVariableID(arg_handler.statusvarname);
        instrument_varID = data->getVariableID(arg_handler.instrumentvarname);
      }

      RelabelingStrategy *relabeling_strategy = new InstrumentalRelabelingStrategy(
              dependent_varID,
              treatment_varID,
              instrument_varID);
      SplittingRule *splitting_rule = new RegressionSplittingRule(data);
      forest = new ForestInstrumental(relabeling_strategy,
                                      splitting_rule,
                                      arg_handler.instrumentvarname);
      break;
    }}

    // Verbose output to logfile if non-verbose mode
    std::ostream* verbose_out;
    if (arg_handler.verbose) {
      verbose_out = &std::cout;
    } else {
      std::ofstream* logfile = new std::ofstream();
      logfile->open("ranger.log");
      if (!logfile->good()) {
        throw std::runtime_error("Could not write to logfile.");
      }
      verbose_out = logfile;
    }

    // Call Ranger
    *verbose_out << "Starting Ranger." << std::endl;
    forest->initCpp(arg_handler.depvarname, arg_handler.memmode, arg_handler.mtry,
        arg_handler.ntree, verbose_out, arg_handler.seed, arg_handler.nthreads,
        arg_handler.predict, arg_handler.targetpartitionsize, arg_handler.splitweights,
        arg_handler.alwayssplitvars, arg_handler.statusvarname, arg_handler.replace,
        arg_handler.savemem,  arg_handler.caseweights, arg_handler.fraction, data);

    forest->run(true);
    if (arg_handler.write) {
      forest->saveToFile();
    }
    forest->writeOutput();
    *verbose_out << "Finished Ranger." << std::endl;

    delete forest;
  } catch (std::exception& e) {
    std::cerr << "Error: " << e.what() << " Ranger will EXIT now." << std::endl;
    delete forest;
    return -1;
  }

  return 0;
}
