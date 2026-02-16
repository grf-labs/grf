/*-------------------------------------------------------------------------------
  Copyright (c) 2024 GRF Contributors.

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

#ifndef GRF_PROGRESSBAR_H_
#define GRF_PROGRESSBAR_H_

#include <atomic>
#include <mutex>
#include <string>

#include "tqdm/tqdm.hpp"
#include "RuntimeContext.h"

namespace grf {

/**
 * Thread-safe wrapper around the tqdm progress bar.
 *
 */
class ProgressBar {
  public:
    ProgressBar(int total,
                const std::string& prefix = "");
    void update(); // called by main thread
    void increment(int n); // called by worked threads
    void finish();

  private:
    int total;
    bool enabled;

    std::atomic<int> done {0};
    std::mutex mtx;
    int last_reported {-1};

    tq::progress_bar pb;
  };

} // namespace grf

#endif // GRF_PROGRESSBAR_H_
