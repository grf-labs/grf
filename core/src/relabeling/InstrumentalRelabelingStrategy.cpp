/*-------------------------------------------------------------------------------
  This file is part of generalized-random-forest (grf).

  generalized-random-forest is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  generalized-random-forest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with generalized-random-forest. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#include "commons/utility.h"
#include "relabeling/InstrumentalRelabelingStrategy.h"

InstrumentalRelabelingStrategy::InstrumentalRelabelingStrategy():
  split_regularization(0) {}

InstrumentalRelabelingStrategy::InstrumentalRelabelingStrategy(double split_regularization):
  split_regularization(split_regularization) {}

std::unordered_map<size_t, double> InstrumentalRelabelingStrategy::relabel(
    const std::vector<size_t>& sampleIDs,
    const Observations& observations) {

  // Prepare the relevant averages.
  size_t num_samples = sampleIDs.size();

  double total_outcome = 0.0;
  double total_treatment = 0.0;
  double total_instrument = 0.0;

  for (size_t sampleID : sampleIDs) {
    total_outcome += observations.get(Observations::OUTCOME, sampleID);
    total_treatment += observations.get(Observations::TREATMENT, sampleID);
    total_instrument += observations.get(Observations::INSTRUMENT, sampleID);
  }

  double average_outcome = total_outcome / num_samples;
  double average_treatment = total_treatment / num_samples;
  double average_instrument = total_instrument / num_samples;
  double average_regularized_instrument = (1 - split_regularization) * average_instrument
    + split_regularization * average_treatment;

  // Calculate the treatment effect.
  double numerator = 0.0;
  double denominator = 0.0;

  for (size_t sampleID : sampleIDs) {
    double outcome = observations.get(Observations::OUTCOME, sampleID);
    double treatment = observations.get(Observations::TREATMENT, sampleID);
    double instrument = observations.get(Observations::INSTRUMENT, sampleID);
    double regularized_instrument = (1 - split_regularization) * instrument + split_regularization * treatment;

    numerator += (regularized_instrument - average_regularized_instrument) * (outcome - average_outcome);
    denominator += (regularized_instrument - average_regularized_instrument) * (treatment - average_treatment);
  }

  if (equal_doubles(denominator, 0.0, 1.0e-10)) {
    return std::unordered_map<size_t, double>(); // Signals that we should not perform a split.
  }

  double local_average_treatment_effect = numerator / denominator;

  // Create the new outcomes.
  std::unordered_map<size_t, double> relabeled_outcomes;

  for (size_t sampleID : sampleIDs) {
    double response = observations.get(Observations::OUTCOME, sampleID);
    double treatment = observations.get(Observations::TREATMENT, sampleID);
    double instrument = observations.get(Observations::INSTRUMENT, sampleID);
    double regularized_instrument = (1 - split_regularization) * instrument + split_regularization * treatment;

    double residual = (response - average_outcome) - local_average_treatment_effect * (treatment - average_treatment);
    relabeled_outcomes[sampleID] = (regularized_instrument - average_regularized_instrument) * residual;
  }
  return relabeled_outcomes;
}
