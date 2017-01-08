#include <map>
#include <unordered_set>
#include <fstream>
#include <random>

#include "catch.hpp"
#include "ForestSerializer.h"
#include "TestUtilities.h"

TEST_CASE("observations serialize and deserialize correctly", "[observationsSerialization]") {
  std::vector<double> outcomes = {-9.99984, -7.36924, 5.11211, -0.826997, 0.655345,
                                  -5.62082, -9.05911, 3.57729, 3.58593, 8.69386};
  std::vector<double> treatment = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::vector<double> instrument = {0, 0, 1, 1, 1, 0, 1, 0, 1, 0};
  Observations original_observations = TestUtilities::create_observations(outcomes, treatment, instrument);

  ObservationsSerializer observations_serializer;
  std::stringstream stream;
  observations_serializer.serialize(stream, original_observations);
  Observations observations = observations_serializer.deserialize(stream);

  REQUIRE(observations.get_num_samples() == original_observations.get_num_samples());

  auto original_observations_by_type = original_observations.get_observations_by_type();
  auto observations_by_type = observations.get_observations_by_type();

  REQUIRE(observations_by_type.size() == original_observations_by_type.size());
  for (auto it = original_observations_by_type.begin(); it != original_observations_by_type.end(); it++) {
    std::string type = it->first;
    REQUIRE(observations_by_type.find(type) != observations_by_type.end());
    REQUIRE(observations_by_type[type].size() == original_observations_by_type[type].size());
  }
}

TEST_CASE("trees serialize and deserialize correctly", "[treeSerialization]") {
  std::shared_ptr<Tree> original_tree(new Tree(
      {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}},
      {{0, 1}, {5, 1}, {4, 6, 9}},
      {10, 20, 30, 40},
      {0.5, 0.75, 0.9, 1.1, 1.2},
      {3, 4, 5, 9, 10, 11},
      {0, 0}));

  TreeSerializer tree_serializer;
  std::stringstream stream;
  tree_serializer.serialize(stream, original_tree);
  std::shared_ptr<Tree> tree = tree_serializer.deserialize(stream);

  REQUIRE(tree->get_sampleIDs().size() == original_tree->get_sampleIDs().size());
  REQUIRE(tree->getChildNodeIDs().size() == original_tree->getChildNodeIDs().size());
  REQUIRE(tree->getSplitVarIDs().size() == original_tree->getSplitVarIDs().size());
  REQUIRE(tree->getSplitValues().size() == original_tree->getSplitValues().size());
  REQUIRE(tree->get_inbag_counts().size() == original_tree->get_inbag_counts().size());
  REQUIRE(tree->getOobSampleIDs().size() == original_tree->getOobSampleIDs().size());
}

TEST_CASE("forests serialize and deserialize correctly", "[forestSerialization]") {
  std::vector<std::shared_ptr<Tree>> trees;
  trees.push_back(std::shared_ptr<Tree>(new Tree({{0}}, {{0}}, {0}, {0}, {0}, {0})));
  trees.push_back(std::shared_ptr<Tree>(new Tree({{1}}, {{1}}, {1}, {1}, {1}, {1})));

  std::vector<double> outcomes = {-9.99984, -7.36924, 5.11211, -0.826997, 0.655345,
                                  -5.62082, -9.05911, 3.57729, 3.58593, 8.69386};
  Observations observations = TestUtilities::create_observations(outcomes);

  ForestSerializer forest_serializer;
  std::stringstream stream;
  Forest original_forest(trees, observations);

  forest_serializer.serialize(stream, original_forest);
  Forest forest = forest_serializer.deserialize(stream);

  REQUIRE(forest.get_trees().size() == original_forest.get_trees().size());
  REQUIRE(forest.get_observations().get_num_samples() == original_forest.get_observations().get_num_samples());
}
