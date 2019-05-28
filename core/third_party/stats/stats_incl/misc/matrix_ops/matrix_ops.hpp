/*################################################################################
  ##
  ##   Copyright (C) 2011-2019 Keith O'Hara
  ##
  ##   This file is part of the StatsLib C++ library.
  ##
  ##   Licensed under the Apache License, Version 2.0 (the "License");
  ##   you may not use this file except in compliance with the License.
  ##   You may obtain a copy of the License at
  ##
  ##       http://www.apache.org/licenses/LICENSE-2.0
  ##
  ##   Unless required by applicable law or agreed to in writing, software
  ##   distributed under the License is distributed on an "AS IS" BASIS,
  ##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  ##   See the License for the specific language governing permissions and
  ##   limitations under the License.
  ##
  ################################################################################*/

/*
 * for internal use only; used to switch between the different matrix libraries
 */

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES

namespace mat_ops
{    
    #include "n_cols.hpp"
    #include "n_rows.hpp"
    #include "n_elem.hpp"

    #include "get_mem_ptr.hpp"

    #include "accu.hpp"
    #include "chol.hpp"
    #include "cumsum.hpp"
    #include "det.hpp"
    #include "exp.hpp"
    #include "eye.hpp"
    #include "fill.hpp"
    #include "get_row.hpp"
    #include "inv.hpp"
    #include "log.hpp"
    #include "log_det.hpp"
    #include "mean.hpp"
    #include "quad_form.hpp"
    #include "repmat.hpp"
    #include "resize.hpp"
    #include "solve.hpp"
    #include "spacing.hpp"
    #include "sum_absdiff.hpp"
    #include "trace.hpp"
    #include "trans.hpp"
    #include "var.hpp"
    #include "zeros.hpp"

    #include "cerr.hpp"
    #include "cout.hpp"
}
#endif
