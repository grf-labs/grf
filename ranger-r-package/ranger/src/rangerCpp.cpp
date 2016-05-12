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

#include <Rcpp.h>
#include <vector>
#include <sstream>

#include "globals.h"
#include "Forest.h"
#include "ForestClassification.h"
#include "ForestRegression.h"
#include "ForestSurvival.h"
#include "ForestProbability.h"
#include "Data.h"
#include "DataChar.h"
#include "DataDouble.h"
#include "DataFloat.h"

// [[Rcpp::export]]
Rcpp::List rangerCpp(uint treetype, std::string dependent_variable_name,
    Rcpp::NumericMatrix input_data, std::vector<std::string> variable_names, uint mtry, uint num_trees, bool verbose,
    uint seed, uint num_threads, bool write_forest, uint importance_mode_r, uint min_node_size,
    std::vector<std::vector<double>>& split_select_weights, bool use_split_select_weights,
    std::vector<std::string>& always_split_variable_names, bool use_always_split_variable_names,
    std::string status_variable_name, bool prediction_mode, Rcpp::List loaded_forest, Rcpp::RawMatrix sparse_data,
    bool sample_with_replacement, bool probability, std::vector<std::string>& unordered_variable_names,
    bool use_unordered_variable_names, bool save_memory, uint splitrule_r, 
    std::vector<double>& case_weights, bool use_case_weights, bool predict_all, 
    bool keep_inbag, double sample_fraction, double alpha, double minprop) {

  Rcpp::List result;
  Forest* forest = 0;
  Data* data = 0;
  try {

    // Empty split select weights and always split variables if not used
    if (!use_split_select_weights) {
      split_select_weights.clear();
    }
    if (!use_always_split_variable_names) {
      always_split_variable_names.clear();
    }
    if (!use_unordered_variable_names) {
      unordered_variable_names.clear();
    }
    if (!use_case_weights) {
      case_weights.clear();
    }

    std::ostream* verbose_out;
    if (verbose) {
      verbose_out = &Rcpp::Rcout;
    } else {
      verbose_out = new std::stringstream;
    }

    size_t num_rows = input_data.nrow();
    size_t num_cols = input_data.ncol();

    // Initialize data with double memmode
    data = new DataDouble(input_data.begin(), variable_names, num_rows, num_cols);

    // If there is sparse data, add it
    if (sparse_data.nrow() > 1) {
      data->addSparseData(sparse_data.begin(), sparse_data.ncol());
    }

    switch (treetype) {
    case TREE_CLASSIFICATION:
      if (probability) {
        forest = new ForestProbability;
      } else {
        forest = new ForestClassification;
      }
      break;
    case TREE_REGRESSION:
      forest = new ForestRegression;
      break;
    case TREE_SURVIVAL:
      forest = new ForestSurvival;
      break;
    case TREE_PROBABILITY:
      forest = new ForestProbability;
      break;
    }

    ImportanceMode importance_mode = (ImportanceMode) importance_mode_r;
    SplitRule splitrule = (SplitRule) splitrule_r;

    // Init Ranger
    forest->initR(dependent_variable_name, data, mtry, num_trees, verbose_out, seed, num_threads,
        importance_mode, min_node_size, split_select_weights, always_split_variable_names, status_variable_name,
        prediction_mode, sample_with_replacement, unordered_variable_names, save_memory, splitrule, case_weights, 
        predict_all, keep_inbag, sample_fraction, alpha, minprop);

    // Load forest object if in prediction mode
    if (prediction_mode) {
      size_t dependent_varID = loaded_forest["dependent.varID"];
      //size_t num_trees = loaded_forest["num.trees"];
      std::vector<std::vector<std::vector<size_t>> > child_nodeIDs = loaded_forest["child.nodeIDs"];
      std::vector<std::vector<size_t>> split_varIDs = loaded_forest["split.varIDs"];
      std::vector<std::vector<double>> split_values = loaded_forest["split.values"];
      std::vector<bool> is_ordered = loaded_forest["is.ordered"];

      if (treetype == TREE_CLASSIFICATION) {
        std::vector<double> class_values = loaded_forest["class.values"];
        ((ForestClassification*) forest)->loadForest(dependent_varID, num_trees, child_nodeIDs, split_varIDs,
            split_values, class_values, is_ordered);
      } else if (treetype == TREE_REGRESSION) {
        ((ForestRegression*) forest)->loadForest(dependent_varID, num_trees, child_nodeIDs, split_varIDs, split_values,
            is_ordered);
      } else if (treetype == TREE_SURVIVAL) {
        size_t status_varID = loaded_forest["status.varID"];
        std::vector<std::vector<std::vector<double>> > chf = loaded_forest["chf"];
        std::vector<double> unique_timepoints = loaded_forest["unique.death.times"];
        ((ForestSurvival*) forest)->loadForest(dependent_varID, num_trees, child_nodeIDs, split_varIDs, split_values,
            status_varID, chf, unique_timepoints, is_ordered);
      } else if (treetype == TREE_PROBABILITY) {
        std::vector<double> class_values = loaded_forest["class.values"];
        std::vector<std::vector<std::vector<double>>>terminal_class_counts =
        loaded_forest["terminal.class.counts"];
        ((ForestProbability*) forest)->loadForest(dependent_varID, num_trees, child_nodeIDs, split_varIDs, split_values,
            class_values, terminal_class_counts, is_ordered);
      }
    }

    // Run Ranger
    forest->run(false);

    if (use_split_select_weights && importance_mode != IMP_NONE) {
      *verbose_out
          << "Warning: Split select weights used. Variable importance measures are only comparable for variables with equal weights."
          << std::endl;
    }

    // Return output
    result.push_back(forest->getNumTrees(), "num.trees");
    result.push_back(forest->getNumIndependentVariables(), "num.independent.variables");
    result.push_back(forest->getPredictions(), "predictions");
    if (treetype == TREE_SURVIVAL) {
      ForestSurvival* temp = (ForestSurvival*) forest;
      result.push_back(temp->getUniqueTimepoints(), "unique.death.times");
    }
    if (!verbose) {
      std::stringstream temp;
      temp << verbose_out->rdbuf();
      result.push_back(temp.str(), "log");
    }
    if (!prediction_mode) {
      result.push_back(forest->getMtry(), "mtry");
      result.push_back(forest->getMinNodeSize(), "min.node.size");
      result.push_back(forest->getVariableImportance(), "variable.importance");
      result.push_back(forest->getOverallPredictionError(), "prediction.error");
    }
    
    if (keep_inbag) {
      result.push_back(forest->getInbagCounts(), "inbag.counts");
    }

    // Save forest if needed
    if (write_forest) {
      Rcpp::List forest_object;
      forest_object.push_back(forest->getDependentVarId(), "dependent.varID");
      forest_object.push_back(forest->getNumTrees(), "num.trees");
      forest_object.push_back(forest->getChildNodeIDs(), "child.nodeIDs");
      forest_object.push_back(forest->getSplitVarIDs(), "split.varIDs");
      forest_object.push_back(forest->getSplitValues(), "split.values");
      forest_object.push_back(forest->getIsOrderedVariable(), "is.ordered");

      if (treetype == TREE_CLASSIFICATION) {
        ForestClassification* temp = (ForestClassification*) forest;
        forest_object.push_back(temp->getClassValues(), "class.values");
      } else if (treetype == TREE_PROBABILITY) {
        ForestProbability* temp = (ForestProbability*) forest;
        forest_object.push_back(temp->getClassValues(), "class.values");
        forest_object.push_back(temp->getTerminalClassCounts(), "terminal.class.counts");
      } else if (treetype == TREE_SURVIVAL) {
        ForestSurvival* temp = (ForestSurvival*) forest;
        forest_object.push_back(temp->getStatusVarId(), "status.varID");
        forest_object.push_back(temp->getChf(), "chf");
        forest_object.push_back(temp->getUniqueTimepoints(), "unique.death.times");
      }
      result.push_back(forest_object, "forest");
    }

    delete forest;
    delete data;
  } catch (std::exception& e) {
    if (strcmp(e.what(), "User interrupt.") != 0) {
      Rcpp::Rcerr << "Error: " << e.what() << " Ranger will EXIT now." << std::endl;
    }
    delete forest;
    delete data;
    return result;
  }

  return result;
}
