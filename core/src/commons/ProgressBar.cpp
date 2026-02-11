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

#include "commons/ProgressBar.h"

namespace grf {

ProgressBar::ProgressBar(int total,
                         std::ostream* progress_bar_output) :
    total(total) {
  if (progress_bar_output == nullptr) {
    pb.set_display(false);
  } else {
    pb.set_ostream(*progress_bar_output);
    pb.set_display(true);
  }
}

void ProgressBar::increment(int n) {
  int v = done.fetch_add(n, std::memory_order_relaxed) + n;
  if (try_lock.try_lock()) {
    pb.update(v, total);
    try_lock.unlock();
  }
}

} // namespace grf
