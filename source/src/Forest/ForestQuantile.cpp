#include <algorithm>
#include <stdexcept>
#include <string>

#include "utility.h"
#include "ForestQuantile.h"
#include "TreeQuantile.h"

ForestQuantile::ForestQuantile(std::vector<double>* quantiles):
    quantiles(quantiles), quantile_predictions(0) {}

ForestQuantile::~ForestQuantile() {}

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
}

void ForestQuantile::growInternal() {
  trees.reserve(num_trees);
  for (size_t i = 0; i < num_trees; ++i) {
    trees.push_back(new TreeQuantile(quantiles));
  }
}

void ForestQuantile::predictInternal() {
  size_t num_prediction_samples = data->getNumRows();
  quantile_predictions.reserve(num_prediction_samples);

  // For all samples get tree predictions
  for (size_t sampleID = 0; sampleID < num_prediction_samples; ++sampleID) {
    if (predict_all) {
      throw std::runtime_error("Quantile forests do not support 'predict_all'.");
    }

    std::unordered_map<size_t, double> weights_by_sampleID;
    for (size_t tree_idx = 0; tree_idx < trees.size(); ++tree_idx) {
      addSampleWeights(sampleID, tree_idx, weights_by_sampleID);
    }
    normalizeSampleWeights(weights_by_sampleID);

    std::vector<double> quantile_cutoffs = calculateQuantileCutoffs(weights_by_sampleID);
    quantile_predictions.push_back(quantile_cutoffs);
  }
}

void ForestQuantile::addSampleWeights(size_t test_sampleID,
                                      size_t tree_idx,
                                      std::unordered_map<size_t, double> &weights_by_sampleID) {
  TreeQuantile *tree = (TreeQuantile *) trees[tree_idx];
  std::vector<size_t> sampleIDs = tree->get_neighboring_samples(test_sampleID);
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

std::vector<double> ForestQuantile::calculateQuantileCutoffs(std::unordered_map<size_t, double> &weights_by_sampleID) {
  std::vector<std::pair<size_t, double>> sampleIDs_and_values;
  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    size_t sampleID = it->first;
    sampleIDs_and_values.push_back(std::pair<size_t, double>(
        sampleID, data->get(sampleID, dependent_varID)));
  }

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
    for (auto it = tree_range.first; it != tree_range.second; ++it) {
      addSampleWeights(sampleID, it->second, weights_by_sampleID);
    }
    normalizeSampleWeights(weights_by_sampleID);

    std::vector<double> quantile_cutoffs = calculateQuantileCutoffs(weights_by_sampleID);
    quantile_predictions.push_back(quantile_cutoffs);
  }
}

void ForestQuantile::writeOutputInternal() {
  *verbose_out << "Tree type:                         " << "Quantile" << std::endl;
}

void ForestQuantile::writeConfusionFile() {

// Open confusion file for writing
  std::string filename = output_prefix + ".confusion";
  std::ofstream outfile;
  outfile.open(filename, std::ios::out);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to confusion file: " + filename + ".");
  }

  *verbose_out << "writing quantile predictions!!!" << std::endl;
  for (auto &quantile_prediction : quantile_predictions) {
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
  std::string filename = output_prefix + ".prediction";
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
  throw std::runtime_error("Saving and loading quantile forests is not yet implemented.");
}

void ForestQuantile::loadFromFileInternal(std::ifstream& infile) {
  throw std::runtime_error("Saving and loading quantile forests is not yet implemented.");
}