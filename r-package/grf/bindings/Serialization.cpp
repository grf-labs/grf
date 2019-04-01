/*-------------------------------------------------------------------------------
  This file is part of generalized random forest (grf).

  grf is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  grf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with grf. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#include <Rcpp.h>
#include <sstream>
#include <vector>

#include "commons/globals.h"
#include "forest/Forest.h"
#include "serialization/ForestSerializer.h"

// [[Rcpp::export]]
Rcpp::RawVector serialize_forest(SEXP xp) {
  Rcpp::XPtr<Forest> forest(xp);
  ForestSerializer forest_serializer;
  std::stringstream stream;
  forest_serializer.serialize(stream, *forest);

  std::string contents = stream.str();

  Rcpp::RawVector result(contents.size());
  std::copy(contents.begin(), contents.end(), result.begin());
  return result;
}

// [[Rcpp::export]]
SEXP deserialize_forest(Rcpp::RawVector input) {
  ForestSerializer forest_serializer;

  std::string contents(input.begin(), input.end());

  std::stringstream stream(contents);
  return Rcpp::XPtr<Forest>(new Forest(forest_serializer.deserialize(stream)));
}
