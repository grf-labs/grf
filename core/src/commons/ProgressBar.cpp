
#include "ProgressBar.h"

ProgressBar::ProgressBar(std::string op_str, size_t max_prog, bool v) {
    operation = std::move(op_str);
    max_progress = max_prog;
    verbose = v;
}

