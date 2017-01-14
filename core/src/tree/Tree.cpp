
#include <iterator>
#include "BootstrapSampler.h"

#include "Tree.h"
#include "utility.h"

Tree::Tree(const std::vector<std::vector<size_t>> &child_nodeIDs,
           const std::vector<std::vector<size_t>>& terminal_nodeIDs,
           const std::vector<size_t>& split_varIDs,
           const std::vector<double>& split_values,
           const std::vector<size_t>& oob_sampleIDs,
           const std::vector<size_t>& inbag_counts) :
    child_nodeIDs(child_nodeIDs),
    terminal_nodeIDs(terminal_nodeIDs),
    split_varIDs(split_varIDs),
    split_values(split_values),
    oob_sampleIDs(oob_sampleIDs),
    inbag_counts(inbag_counts) {}

Tree::~Tree() {}


