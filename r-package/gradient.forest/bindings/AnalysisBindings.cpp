/*-------------------------------------------------------------------------------
  This file is part of generalized-random-forest.

  gradient-forest is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  gradient-forest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with gradient-forest. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

// [[Rcpp::export]]
Rcpp::NumericMatrix compute_split_frequencies(Rcpp::List forest,
                                             size_t max_depth) {
  Forest deserialized_forest = RcppUtilities::deserialize_forest(
      forest[RcppUtilities::SERIALIZED_FOREST_KEY]);

  VariableImportanceComputer computer;
  std::vector<std::vector<size_t>> variable_frequencies = computer.compute(forest, max_depth);

  size_t num_variables = forest->get_num_variables;
  Rcpp::NumericMatrix result(num_variables, max_depth);
  for (size_t i = 0; i < num_variables; i++) {
    const std::vector<double>& frequencies = variable_frequencies.at(i);
    for (size_t j = 0; j < frequencies.size(); j++) {
      double frequency = frequencies[j];
      result(i, j) = frequency;
    }
  }
  return result;
}
