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
                         const std::string& prefix) :
    total(total) {
  std::ostream* out = grf::runtime_context.verbose_stream;
  if (out == nullptr) {
    enabled = false;
  } else {
    pb.set_display(true);
    pb.set_ostream(*out);
    pb.set_prefix(prefix);
    pb.set_bar_symbol("\033[38;5;65m\u2588\033[0m"); // grf forest-greenish color.
    enabled = true;
  }
}

void ProgressBar::increment(int n) {
  if (!enabled) return;

  int current = done.fetch_add(n, std::memory_order_relaxed) + n;
  if (mtx.try_lock()) {
    pb.update(current, total);
    last_reported = current;
    mtx.unlock();
  }
}

// Only do a final update if the progress bar hasn't already been updated to the total.
void ProgressBar::finish() {
  if (!enabled) return;

  std::lock_guard<std::mutex> lock(mtx); // final call happens outside the multi-threaded context, but we keep the lock just for safety.
  if (last_reported < total) {
    pb.update(total, total);
  }
}

} // namespace grf
