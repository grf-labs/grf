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
#include <cmath>

#include "commons/utility.h"
#include "ObjectiveBayesDebiaser.h"

namespace grf {

double ObjectiveBayesDebiaser::debias(double var_between,
                                      double group_noise,
                                      double num_good_groups) const {
  
  // Let S denote the true between-groups variance, and assume that
  // group_noise is measured exactly; our method-of-moments estimate is
  // then \hat{S} = var_between - group_noise. Now, if we take
  // num_good_groups * var_between to be chi-squared with scale
  // S + group_noise and num_good_groups degrees of freedom, and assume
  // that group_noise >> S, we find that var[initial_estimate] is roughly
  // var_between^2 * 2 / num_good_groups; moreover, the distribution of
  // \hat{S} - S is roughly Gaussian with this variance. Our estimation strategy
  // relies on this fact, and puts a uniform prior on S for the interval [0, infty).
  // This debiasing does nothing when \hat{S} >> var_between * sqrt(2 / num_good_groups),
  // but keeps \hat{S} from going negative.
  
  double initial_estimate = var_between - group_noise;
  double initial_se = std::max(var_between, group_noise) * std::sqrt(2.0 / num_good_groups);
  
  double ratio = initial_estimate / initial_se;
  
  // corresponds to \int_(-r)^infty x * phi(x) dx, for the standard Gaussian density phi(x)
  double numerator = std::exp(- ratio * ratio / 2) * ONE_over_SQRT_TWO_PI;

  // corresponds to int_(-r)^infty phi(x) dx
  double denominator =  0.5 * std::erfc(-ratio * ONE_over_SQRT_TWO);

  // this is the o-Bayes estimate of the error of the initial estimate
  double bayes_correction = initial_se * numerator / denominator;
  return initial_estimate + bayes_correction;
}

} // namespace grf
