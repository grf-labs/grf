/*-------------------------------------------------------------------------------
  This file is part of gradient-forest.

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

#ifndef GRADIENTFOREST_OBJECTBAYESDEBIASER_H
#define GRADIENTFOREST_OBJECTBAYESDEBIASER_H


class ObjectiveBayesDebiaser {
public:
  double debias(double within_noise,
                double num_good_groups,
                double var_between);
};

#endif //GRADIENTFOREST_OBJECTBAYESDEBIASER_H
