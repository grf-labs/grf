#include <algorithm>
#include <stdexcept>
#include <string>

#include "utility.h"
#include "ForestQuantile.h"

ForestQuantile::ForestQuantile(std::vector<double>* quantiles):
    quantiles(quantiles), original_responses(0) {}

ForestQuantile::~ForestQuantile() {}

void ForestQuantile::loadForest(size_t dependent_varID, size_t num_trees,
                                std::vector<std::vector<std::vector<size_t>>> &forest_child_nodeIDs,
                                std::vector<std::vector<size_t>> &forest_split_varIDs,
                                std::vector<std::vector<double>> &forest_split_values,
                                std::vector<double>* quantiles,
                                std::vector<std::vector<std::vector<size_t>>> sampleIDs,
                                std::vector<double>* originalResponses) {

  this->dependent_varID = dependent_varID;
  this->num_trees = num_trees;
  this->original_responses = new std::vector<double>(*originalResponses);

  // Create trees
  trees.reserve(num_trees);
  for (size_t i = 0; i < num_trees; ++i) {
    Tree* tree = new TreeQuantile(forest_child_nodeIDs[i], forest_split_varIDs[i], forest_split_values[i],
                                  quantiles, sampleIDs[i]);
    trees.push_back(tree);
  }

  // Create thread ranges
  equalSplit(thread_ranges, 0, num_trees - 1, num_threads);
}

void ForestQuantile::initInternal(std::string status_variable_name) {

  // If mtry not set, use number of independent variables / 3.
  if (mtry == 0) {
    unsigned long temp = sqrt((double) (num_variables - 1));
    mtry = std::max((unsigned long) 1, temp);
  }

  // Set minimal node size
  if (min_node_size == 0) {
    min_node_size = DEFAULT_MIN_NODE_SIZE_REGRESSION;
  }

  // Sort data if memory saving mode
  if (!memory_saving_splitting) {
    data->sort();
  }

  if (!prediction_mode) {
    original_responses = new std::vector<double>();
    for (size_t sampleID = 0; sampleID < data->getNumRows(); ++sampleID) {
      original_responses->push_back(data->get(sampleID, dependent_varID));
    }
  }
}

void ForestQuantile::growInternal() {
  trees.reserve(num_trees);
  for (size_t i = 0; i < num_trees; ++i) {
    trees.push_back(new TreeQuantile(quantiles));
  }
}

void ForestQuantile::predictInternal() {
  if (predict_all) {
    throw std::runtime_error("Quantile forests do not support 'predict_all'.");
  }
  
  for (size_t sampleID = 0; sampleID < data->getNumRows(); ++sampleID) {
    std::unordered_map<size_t, double> weights_by_sampleID;
    for (size_t tree_idx = 0; tree_idx < trees.size(); ++tree_idx) {
      TreeQuantile* tree = (TreeQuantile *) trees[tree_idx];
      addSampleWeights(sampleID, tree, weights_by_sampleID);
    }

    normalizeSampleWeights(weights_by_sampleID);

    std::vector<std::pair<size_t, double>> sampleIDs_and_values;
    for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
      size_t sampleID = it->first;
      sampleIDs_and_values.push_back(std::pair<size_t, double>(
          sampleID, original_responses->at(sampleID)));
    }

    std::vector<double> quantile_cutoffs = calculateQuantileCutoffs(weights_by_sampleID,
                                                                    sampleIDs_and_values);

    predictions.push_back(quantile_cutoffs);
  }
}

void ForestQuantile::addSampleWeights(size_t test_sample_idx,
                                      TreeQuantile* tree,
                                      std::unordered_map<size_t, double> &weights_by_sampleID) {
  std::vector<size_t> sampleIDs = tree->get_neighboring_samples(test_sample_idx);
  double sample_weight = 1.0 / sampleIDs.size();

  for (auto &sampleID : sampleIDs) {
    weights_by_sampleID[sampleID] += sample_weight;
  }
}

void ForestQuantile::normalizeSampleWeights(std::unordered_map<size_t, double>& weights_by_sampleID) {
  double total_weight = 0.0;
  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    total_weight += it->second;
  }

  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    it->second /= total_weight;
  }
}

std::vector<double> ForestQuantile::calculateQuantileCutoffs(std::unordered_map<size_t,double> &weights_by_sampleID,
                                                             std::vector<std::pair<size_t, double>> sampleIDs_and_values) {
  std::sort(sampleIDs_and_values.begin(),
            sampleIDs_and_values.end(),
            [](std::pair<size_t, double> first_pair, std::pair<size_t, double> second_pair) {
              return first_pair.second < second_pair.second;
            });

  std::vector<double> quantile_cutoffs;
  auto quantile_it = quantiles->begin();
  double cumulative_weight = 0.0;

  for (auto it = sampleIDs_and_values.begin(); it != sampleIDs_and_values.end()
                                               && quantile_it != quantiles->end(); ++it) {
    size_t sampleID = it->first;
    double value = it->second;

    cumulative_weight += weights_by_sampleID[sampleID];
    if (cumulative_weight >= *quantile_it) {
      quantile_cutoffs.push_back(value);
      ++quantile_it;
    }
  }

  double last_value = sampleIDs_and_values.back().first;
  for (; quantile_it != quantiles->end(); ++quantile_it) {
    quantile_cutoffs.push_back(last_value);
  }
  return quantile_cutoffs;
}

void ForestQuantile::computePredictionErrorInternal() {
  if (predict_all) {
    throw std::runtime_error("Quantile forests do not support 'predict_all'.");
  }

  std::unordered_multimap<size_t, size_t> trees_by_oob_samples;
  for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
    for (size_t sampleID : trees[tree_idx]->getOobSampleIDs()) {
      trees_by_oob_samples.insert(std::pair<size_t, size_t>(sampleID, tree_idx));
    }
  }

  for (size_t sampleID = 0; sampleID < data->getNumRows(); ++sampleID) {
    std::unordered_map<size_t, double> weights_by_sampleID;

    auto tree_range = trees_by_oob_samples.equal_range(sampleID);
    if (tree_range.first == tree_range.second) {
      predictions.push_back(std::vector<double>(quantiles->size(), 0.0));
      continue;
    }

    for (auto it = tree_range.first; it != tree_range.second; ++it) {
      TreeQuantile* tree = (TreeQuantile*) trees[it->second];

      // hackhackhack
      std::vector<size_t> oob_sampleIDs = tree->getOobSampleIDs();
      size_t sample_idx = (size_t) (std::find(oob_sampleIDs.begin(), oob_sampleIDs.end(), sampleID) - oob_sampleIDs.begin());

      addSampleWeights(sample_idx, tree, weights_by_sampleID);
    }

    normalizeSampleWeights(weights_by_sampleID);

    std::vector<std::pair<size_t, double>> sampleIDs_and_values;
    for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
      size_t sampleID = it->first;
      sampleIDs_and_values.push_back(std::pair<size_t, double>(
          sampleID, data->get(sampleID, dependent_varID)));
    }

    std::vector<double> quantile_cutoffs = calculateQuantileCutoffs(weights_by_sampleID,
                                                                    sampleIDs_and_values);
    predictions.push_back(quantile_cutoffs);
  }
}

void ForestQuantile::writeOutputInternal() {
  *verbose_out << "Tree type:                         " << "Quantile" << std::endl;
}

void ForestQuantile::writeConfusionFile() {

// Open confusion file for writing
  std::string filename = "ranger.confusion";
  std::ofstream outfile;
  outfile.open(filename, std::ios::out);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to confusion file: " + filename + ".");
  }

  for (auto &quantile_prediction : predictions) {
    for (auto it = quantile_prediction.begin(); it != quantile_prediction.end(); ++it) {
      outfile << *it << " ";
    }
    outfile << std::endl;
  }

  outfile.close();
  *verbose_out << "Saved prediction error to file " << filename << "." << std::endl;
}

void ForestQuantile::writePredictionFile() {

// Open prediction file for writing
  std::string filename = "ranger.prediction";
  std::ofstream outfile;
  outfile.open(filename, std::ios::out);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to prediction file: " + filename + ".");
  }

  // Write
  outfile << "Predictions: " << std::endl;
  for (size_t i = 0; i < predictions.size(); ++i) {
    for (size_t j = 0; j < predictions[i].size(); ++j) {
      outfile << predictions[i][j] << " ";
    }
    outfile << std::endl;
  }

  *verbose_out << "Saved predictions to file " << filename << "." << std::endl;
}

void ForestQuantile::saveToFileInternal(std::ofstream& outfile) {

// Write num_variables
  outfile.write((char*) &num_variables, sizeof(num_variables));

// Write treetype
  TreeType treetype = TREE_QUANTILE;
  outfile.write((char*) &treetype, sizeof(treetype));
  saveVector1D(*quantiles, outfile);
}

void ForestQuantile::loadFromFileInternal(std::ifstream& infile) {

// Read number of variables
  size_t num_variables_saved;
  infile.read((char *) &num_variables_saved, sizeof(num_variables_saved));

// Read treetype
  TreeType treetype;
  infile.read((char *) &treetype, sizeof(treetype));
  if (treetype != TREE_QUANTILE) {
    throw std::runtime_error("Wrong treetype. Loaded file is not a regression forest.");
  }

  readVector1D(*quantiles, infile);

  for (size_t i = 0; i < num_trees; ++i) {

    // Read data
    std::vector<std::vector<size_t>> child_nodeIDs;
    readVector2D(child_nodeIDs, infile);
    std::vector<size_t> split_varIDs;
    readVector1D(split_varIDs, infile);
    std::vector<double> split_values;
    readVector1D(split_values, infile);

    // If dependent variable not in test data, change variable IDs accordingly
    if (num_variables_saved > num_variables) {
      for (auto &varID : split_varIDs) {
        if (varID >= dependent_varID) {
          --varID;
        }
      }
    }

    std::vector<std::vector<size_t>> sampleIDs;
    readVector2D(sampleIDs, infile);

    // Create tree
    Tree *tree = new TreeQuantile(child_nodeIDs, split_varIDs, split_values, quantiles, sampleIDs);
    trees.push_back(tree);
  }
}