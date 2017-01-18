#include <vector>
#include <string>
#include "Observations.h"
#include "utility.h"
#include "InstrumentalPredictionStrategy.h"

size_t InstrumentalPredictionStrategy::prediction_length() {
    return 1;
}

Prediction InstrumentalPredictionStrategy::predict(const std::map<std::string, double>& average_prediction_values,
                                                   const std::unordered_map<size_t, double>& weights_by_sampleID,
                                                   const Observations& observations) {
  // Compute the relevant averages.
  double average_instrument = 0.0;
  double average_treatment = 0.0;
  double average_outcome = 0.0;

  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    size_t neighborID = it->first;
    double weight = it->second;

    average_outcome += weight * observations.get(Observations::OUTCOME)[neighborID];
    average_treatment += weight * observations.get(Observations::TREATMENT)[neighborID];
    average_instrument += weight * observations.get(Observations::INSTRUMENT)[neighborID];
  }

  // Finally, calculate the prediction.
  double instrument_effect = 0.0;
  double first_stage_effect = 0.0;
  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    size_t neighborID = it->first;
    double weight = it->second;

    double response = observations.get(Observations::OUTCOME)[neighborID];
    double treatment = observations.get(Observations::TREATMENT)[neighborID];
    double instrument = observations.get(Observations::INSTRUMENT)[neighborID];

    instrument_effect += weight * (instrument - average_instrument) * (response - average_outcome);
    first_stage_effect += weight * (instrument - average_instrument) * (treatment - average_treatment);
  }

  return Prediction({ instrument_effect / first_stage_effect });
}

Prediction InstrumentalPredictionStrategy::predict_with_variance(
    const std::vector<std::vector<size_t>>& leaf_sampleIDs,
    const Observations& observations,
    uint ci_group_size) {

  size_t num_trees = leaf_sampleIDs.size();

  std::vector<double> leaf_Y(num_trees);
  std::vector<double> leaf_W(num_trees);
  std::vector<double> leaf_Z(num_trees);
  std::vector<double> leaf_YZ(num_trees);
  std::vector<double> leaf_WZ(num_trees);

  double total_Y = 0;
  double total_W = 0;
  double total_Z = 0;
  double total_YZ = 0;
  double total_WZ = 0;

  double non_empty_leaves = 0;

  for (int i = 0; i < leaf_sampleIDs.size(); ++i) {

    size_t leaf_size = leaf_sampleIDs[i].size();

    if (leaf_size == 0) {

      leaf_Y[i] = NAN;
      leaf_W[i] = NAN;
      leaf_Z[i] = NAN;
      leaf_YZ[i] = NAN;
      leaf_WZ[i] = NAN;

      continue;

    }

    double sum_Y = 0;
    double sum_W = 0;
    double sum_Z = 0;
    double sum_YZ = 0;
    double sum_WZ = 0;

    for (auto sampleID :leaf_sampleIDs[i]) {
      sum_Y += observations.get(Observations::OUTCOME)[sampleID];
      sum_W += observations.get(Observations::TREATMENT)[sampleID];
      sum_Z += observations.get(Observations::INSTRUMENT)[sampleID];
      sum_YZ += observations.get(Observations::OUTCOME)[sampleID] * observations.get(Observations::INSTRUMENT)[sampleID];
      sum_WZ += observations.get(Observations::TREATMENT)[sampleID] * observations.get(Observations::INSTRUMENT)[sampleID];
    }

    leaf_Y[i] = sum_Y / leaf_size;
    leaf_W[i] = sum_W / leaf_size;
    leaf_Z[i] = sum_Z / leaf_size;
    leaf_YZ[i] = sum_YZ / leaf_size;
    leaf_WZ[i] = sum_WZ / leaf_size;

    total_Y += sum_Y / leaf_size;
    total_W += sum_W / leaf_size;
    total_Z += sum_Z / leaf_size;
    total_YZ += sum_YZ / leaf_size;
    total_WZ += sum_WZ / leaf_size;

    non_empty_leaves++;

  }

  double avg_Y = total_Y / non_empty_leaves;
  double avg_W = total_W / non_empty_leaves;
  double avg_Z = total_Z / non_empty_leaves;
  double avg_YZ = total_YZ / non_empty_leaves;
  double avg_WZ = total_WZ / non_empty_leaves;

  double instrument_effect = avg_YZ - avg_Y * avg_Z;
  double first_stage = avg_WZ - avg_W * avg_Z;
  double treatment_estimate = instrument_effect / first_stage;
  double main_effect = avg_Y - avg_W * treatment_estimate;

  double num_good_groups = 0;
  std::vector<double> psi_mean = {0, 0}; // need to center this over good groups
  std::vector<std::vector<double>> psi_squared = {{0, 0}, {0, 0}};
  std::vector<std::vector<double>> psi_grouped_squared = {{0, 0}, {0, 0}};

  for (int group = 0; group < leaf_sampleIDs.size() / ci_group_size; ++group) {

    bool good_group = true;
    for (int j = 0; j < ci_group_size; ++j) {
      if (leaf_sampleIDs[group * ci_group_size + j].size() == 0) {
        good_group = false;
      }
    }
    if (!good_group) continue;

    num_good_groups++;

    double group_psi_1 = 0;
    double group_psi_2 = 0;

    for (int j = 0; j < ci_group_size; ++j) {

      int i = group * ci_group_size + j;
      double psi_1 = leaf_YZ[i] - leaf_WZ[i] * treatment_estimate - leaf_Z[i] * main_effect;
      double psi_2 = leaf_Y[i] - leaf_W[i] * treatment_estimate - main_effect;

      psi_squared[0][0] += psi_1 * psi_1;
      psi_squared[0][1] += psi_1 * psi_2;
      psi_squared[1][0] += psi_2 * psi_1;
      psi_squared[1][1] += psi_2 * psi_2;

      group_psi_1 += psi_1;
      group_psi_2 += psi_2;
    }

    group_psi_1 /= ci_group_size;
    group_psi_2 /= ci_group_size;

    psi_mean[0] += group_psi_1;
    psi_mean[1] += group_psi_2;

    psi_grouped_squared[0][0] += group_psi_1 * group_psi_1;
    psi_grouped_squared[0][1] += group_psi_1 * group_psi_2;
    psi_grouped_squared[1][0] += group_psi_2 * group_psi_1;
    psi_grouped_squared[1][1] += group_psi_2 * group_psi_2;
  }

  for (int i = 0; i < 2; i++) {
    psi_mean[i] /= num_good_groups;
    for (int j = 0; j < 2; j++) {
      psi_squared[i][j] /= (num_good_groups * ci_group_size);
      psi_grouped_squared[i][j] /= num_good_groups;
    }
  }

  std::vector<std::vector<double>> psi_half_sample_variance = {{0, 0}, {0, 0}};

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      psi_half_sample_variance[i][j] = psi_grouped_squared[i][j]
          - (psi_squared[i][j] - psi_grouped_squared[i][j]) / (ci_group_size - 1)
          - psi_mean[i] * psi_mean[j];
    }
  }

  double matrix_product = psi_half_sample_variance[0][0]
                              - psi_half_sample_variance[0][1] * avg_W
                              - psi_half_sample_variance[1][0] * avg_W
                              + psi_half_sample_variance[1][1] * avg_W * avg_W;
  double variance_estimate = matrix_product / (first_stage * first_stage);

  return Prediction({ treatment_estimate }, { variance_estimate });
}

bool InstrumentalPredictionStrategy::requires_leaf_sampleIDs() {
  return true;
}

PredictionValues InstrumentalPredictionStrategy::precompute_prediction_values(
    const std::vector<std::vector<size_t>> &leaf_sampleIDs,
    const Observations &observations) {
  return PredictionValues();
}
