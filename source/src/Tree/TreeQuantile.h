#ifndef RANGER_TREEQUANTILE_H
#define RANGER_TREEQUANTILE_H

#include "TreeRegression.h"
#include "TreeClassification.h"

class TreeQuantile: public TreeRegression {
public:
  TreeQuantile(std::vector<double>* quantiles);
  bool splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs);
  std::vector<size_t> get_neighboring_samples(size_t sampleID);

private:
  std::vector<uint>* relabelResponses(std::vector<double>* responses);
  TreeClassification* createClassificationTree(std::vector<size_t>& nodeSampleIDs,
                                               std::vector<uint>* relabeledResponses);

  std::vector<double>* quantiles;

  DISALLOW_COPY_AND_ASSIGN(TreeQuantile);
};
#endif //RANGER_TREEQUANTILE_H
