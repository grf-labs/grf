#ifndef RANGER_TreeInstrumental_H
#define RANGER_TreeInstrumental_H

#include "RegressionTreeFactory.h"

class InstrumentalTreeFactory: public TreeFactory {
public:
  InstrumentalTreeFactory(size_t treatment_varID, size_t instrument_varID, std::string instrument_var_name);

  InstrumentalTreeFactory(std::vector<std::vector<size_t>> &child_nodeIDs, std::vector<size_t> &split_varIDs,
                   std::vector<double> &split_values,
                   std::vector<std::vector<size_t>> sampleIDs,
                   size_t treatment_varID, size_t instrument_varID,
                   std::string instrument_var_name);
  std::vector<size_t> get_neighboring_samples(size_t sampleID);

  bool splitNodeInternal(size_t nodeID, std::vector<size_t> &possible_split_varIDs);

private:
  size_t treatment_varID;
  size_t instrument_varID;
  std::string instrument_var_name;

  DISALLOW_COPY_AND_ASSIGN(InstrumentalTreeFactory);
};
#endif //RANGER_TreeInstrumental_H
