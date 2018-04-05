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

#ifndef GRF_INSTRUMENTALRELABELINGSTRATEGY_H
#define GRF_INSTRUMENTALRELABELINGSTRATEGY_H

#include <unordered_map>
#include <vector>
#include "commons/Observations.h"
#include "tree/Tree.h"
#include "relabeling/RelabelingStrategy.h"

class InstrumentalRelabelingStrategy: public RelabelingStrategy {
public:
  InstrumentalRelabelingStrategy();

  InstrumentalRelabelingStrategy(double reduced_form_weight);

  std::unordered_map<size_t, double> relabel(
      const std::vector<size_t>& samples,
      const Observations& observations);

  DISALLOW_COPY_AND_ASSIGN(InstrumentalRelabelingStrategy);

private:
  double reduced_form_weight;
};

#endif //GRF_INSTRUMENTALRELABELINGSTRATEGY_H
