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
#------------------------------------------------------------------------  */

#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H
#include <cassert>
#include <random>
#include <vector>
#include <array>

namespace distributions {

    template<typename T>
    T uniform_draw(const double min_value, const double max_value, std::mt19937_64 &eng) {
        for (;;) {
            auto numerator = static_cast<double>(eng() - (eng.min)());
            auto divisor = static_cast<double>((eng.max)() - (eng.min)());
            assert(divisor > 0);
            assert(numerator >= 0 && numerator <= divisor);
            double result = numerator / divisor * (max_value - min_value) + min_value;
            if (result < max_value) return static_cast<T>(result);
        }
    }

    template<typename T>
    T poisson_draw(const double mean, std::mt19937_64 &eng) {
        using std::floor;
        using std::abs;
        using std::log;

        double _expmean = exp(-mean);

        // Use simple method if mean is small enough
        if (mean < 10) {
            double p = _expmean;
            double x = 0;
            auto u = uniform_draw<double>(0., 1., eng);
            while (u > p) {
                u = u - p;
                ++x;
                p = mean * p / x;
            }
            return x;
        } else {
            // Use more complicated method for larger means

            double smu = sqrt(mean);
            double b = 0.931 + 2.53 * smu;
            double a = -0.059 + 0.02483 * b;
            double inv_alpha = 1.1239 + 1.1328 / (b - 3.4);
            double v_r = 0.9277 - 3.6224 / (b - 2);
            double log_sqrt_2pi = 0.91893853320467267;
            std::array<double, 10> poisson_table = {
                    0.0,
                    0.0,
                    0.69314718055994529,
                    1.7917594692280550,
                    3.1780538303479458,
                    4.7874917427820458,
                    6.5792512120101012,
                    8.5251613610654147,
                    10.604602902745251,
                    12.801827480081469
            };


            while (true) {
                double u;
                auto v = uniform_draw<double>(0., 1., eng);
                if (v <= 0.86 * v_r) {
                    u = v / v_r - 0.43;
                    return static_cast<double>(floor((2 * a / (0.5 - abs(u)) + b) * u + mean + 0.445));
                }

                if (v >= v_r) {
                    u = uniform_draw<double>(0., 1., eng) - 0.5;
                } else {
                    u = v / v_r - 0.93;
                    u = ((u < 0) ? -0.5 : 0.5) - u;
                    v = uniform_draw<double>(0., 1., eng) * v_r;
                }

                double us = 0.5 - abs(u);
                if (us < 0.013 && v > us) {
                    continue;
                }

                double k = floor((2 * a / us + b) * u + mean + 0.445);
                v = v * inv_alpha / (a / (us * us) + b);

                if (k >= 10) {
                    if (log(v * smu) <= (k + 0.5) * log(mean / k)
                                        - mean
                                        - log_sqrt_2pi
                                        + k
                                        - (1 / 12. - (1 / 360. - 1 / (1260. * k * k)) / (k * k)) / k) {
                        return static_cast<T>(k);
                    }
                } else if (k >= 0) {
                    if (log(v) <= k * log(mean) - mean - poisson_table[static_cast<size_t>(k)]) {
                        return static_cast<T>(k);
                    }
                }
            }
        }
    }

}



#endif // DISTRIBUTIONS_H
