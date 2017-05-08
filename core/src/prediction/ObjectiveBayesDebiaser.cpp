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

#include "commons/utility.h"
#include "ObjectiveBayesDebiaser.h"


double ObjectiveBayesDebiaser::debias(double within_noise,
                                      double num_good_groups,
                                      double var_between) {
  // Start by computing chi-squared log-density, taking within_noise as fixed.
  double lx[100];
  double x[100];
  for (size_t iter = 0; iter < 100; ++iter) {
    x[iter] = iter * within_noise / 400.0;
    lx[iter] = std::log(x[iter] + within_noise) * (1.0 - num_good_groups / 2.0)
               - (num_good_groups * var_between) / (2 * (x[iter] + within_noise));
  }

  // Compute maximal log-density, to stabilize call to exp() below.
  double maxlx = lx[1];
  for (size_t iter = 0; iter < 100; ++iter) {
    maxlx = std::max(maxlx, lx[iter]);
  }

  // Now do Bayes rule.
  double numerator = 0;
  double denominator = 0;
  for (size_t iter = 0; iter < 100; ++iter) {
    double fx = std::exp(lx[iter] - maxlx);
    numerator += x[iter] * fx;
    denominator += fx;
  }

  double result = numerator / denominator;

  // Avoid jumpy behavior at the threshold.
  return std::min(result, within_noise / 10);
}