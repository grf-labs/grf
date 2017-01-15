#include "Tree.h"
#include "TreeSerializer.h"
#include "utility.h"
#include "PredictionValuesSerializer.h"

void TreeSerializer::serialize(std::ostream& stream, std::shared_ptr<Tree> tree) {
  saveVector2D(tree->get_child_nodeIDs(), stream);
  saveVector2D(tree->get_leaf_nodeIDs(), stream);
  saveVector1D(tree->get_split_varIDs(), stream);
  saveVector1D(tree->get_split_values(), stream);
  saveVector1D(tree->get_oob_sampleIDs(), stream);
  saveVector1D(tree->get_inbag_counts(), stream);

  PredictionValuesSerializer prediction_values_serializer;
  prediction_values_serializer.serialize(stream, tree->get_prediction_values());
}

std::shared_ptr<Tree> TreeSerializer::deserialize(std::istream& stream) {
  std::vector<std::vector<size_t>> child_nodeIDs;
  readVector2D(child_nodeIDs, stream);

  std::vector<std::vector<size_t>> sampleIDs;
  readVector2D(sampleIDs, stream);

  std::vector<size_t> split_varIDs;
  readVector1D(split_varIDs, stream);

  std::vector<double> split_values;
  readVector1D(split_values, stream);

  std::vector<size_t> oob_sampleIDs;
  readVector1D(oob_sampleIDs, stream);

  std::vector<size_t> inbag_counts;
  readVector1D(inbag_counts, stream);

  PredictionValuesSerializer prediction_values_serializer;
  PredictionValues prediction_values = prediction_values_serializer.deserialize(stream);

  return std::shared_ptr<Tree>(
      new Tree(child_nodeIDs,
               sampleIDs,
               split_varIDs,
               split_values,
               oob_sampleIDs,
               inbag_counts,
               prediction_values));
}
