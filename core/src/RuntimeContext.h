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

#ifndef GRF_RUNTIME_CONTEXT_H_
#define GRF_RUNTIME_CONTEXT_H_

#include <functional>
#include <ostream>
#include <string>

namespace grf {

/**
 * Runtime context for language bindings.
 *
 * - forest_name: Used in progress bar labels
 * - verbose_stream: Output stream for verbose logging (nullptr disables)
 * - interrupt_handler: Called periodically during long operations to check for
 *   user interrupts (e.g., Ctrl-C). Language bindings should set this to their
 *   interrupt check function (R: Rcpp::checkUserInterrupt, Python: PyErr_CheckSignals).
 *   Default is a no-op for standalone C++ usage.
 */
struct RuntimeContext {
  std::string forest_name = "grf";
  std::ostream* verbose_stream = nullptr;
  std::function<void()> interrupt_handler = []() {};
};

extern RuntimeContext runtime_context;

} // namespace grf

#endif // GRF_RUNTIME_CONTEXT_H_
