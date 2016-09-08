#ifndef RANGER_TreeInstrumental_H
#define RANGER_TreeInstrumental_H

#include "TreeRegression.h"
#include "TreeClassification.h"

class TreeInstrumental: public TreeRegression {
public:
  TreeInstrumental(size_t treatment_varID, size_t instrument_varID, std::string instrument_var_name);

  TreeInstrumental(std::vector<std::vector<size_t>> &child_nodeIDs, std::vector<size_t> &split_varIDs,
                   std::vector<double> &split_values, std::vector<bool> *is_ordered_variable,
                   std::vector<std::vector<size_t>> sampleIDs,
                   size_t treatment_varID, size_t instrument_varID,
                   std::string instrument_var_name);
  std::vector<size_t> get_neighboring_samples(size_t sampleID);

  bool splitNodeInternal(size_t nodeID, std::vector<size_t> &possible_split_varIDs);

private:
  std::unordered_map<size_t, double> relabelResponses(std::vector<size_t>& nodeSampleIDs);
  bool equalDoubles(double first, double second);
  int sgn(double val);

  void appendToFileInternal(std::ofstream& file);

  size_t treatment_varID;
  size_t instrument_varID;
  std::string instrument_var_name;
  std::uniform_int_distribution<uint> udist;

  DISALLOW_COPY_AND_ASSIGN(TreeInstrumental);
};
#endif //RANGER_TreeInstrumental_H
