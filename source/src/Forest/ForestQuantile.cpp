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
  for (size_t sample_idx = 0; sample_idx < num_prediction_samples; ++sample_idx) {
    if (predict_all) {
      throw std::runtime_error("Quantile forests do not support 'predict_all'.");
    }

    std::unordered_map<size_t, double> weights_by_sampleID = calculateSampleWeights(sample_idx);
    std::vector<double> quantile_cutoffs = calculateQuantileCutoffs(weights_by_sampleID);
    quantile_predictions.push_back(quantile_cutoffs);
  }
}

std::unordered_map<size_t, double> ForestQuantile::calculateSampleWeights(size_t sample_idx) {
  std::unordered_map<size_t, double> weights_by_sampleID;

  // Weight each training sample by how often it appears in the same terminal
  // node as the test sample.
  double total_weight = 0.0;
  for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
    TreeQuantile* tree = (TreeQuantile*) trees[tree_idx];

    std::vector<size_t> sampleIDs = tree->get_neighboring_samples(sample_idx);
    double sample_weight = 1.0 / sampleIDs.size();

    std::cout << "sample_weight " << sample_weight << std::endl;
    for (auto& sampleID : sampleIDs) {
      weights_by_sampleID[sampleID] += sample_weight;
      total_weight += sample_weight;
    }
  }

  // Normalize the weights so that they sum to one.
  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    it->second /= total_weight;
  }

  return weights_by_sampleID;
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
  for (int sample_idx : {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}) {
    std::unordered_map<size_t, double> weights_by_sampleID = calculateSampleWeights(sample_idx);
    std::vector<double> quantile_cutoffs = calculateQuantileCutoffs(weights_by_sampleID);
    quantile_predictions.push_back(quantile_cutoffs);
  }

  std::cout << "quantile predictions!!!" << std::endl;
  for (auto &quantile_prediction : quantile_predictions) {
    for (auto it = quantile_prediction.begin(); it != quantile_prediction.end(); ++it) {
      std::cout << *it << " ";
    }
    std::cout << std::endl;
  }
}

void ForestQuantile::writeOutputInternal() {
  *verbose_out << "Tree type:                         " << "Regression" << std::endl;
}

void ForestQuantile::writeConfusionFile() {

// Open confusion file for writing
  std::string filename = output_prefix + ".confusion";
  std::ofstream outfile;
  outfile.open(filename, std::ios::out);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to confusion file: " + filename + ".");
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