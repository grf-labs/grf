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
 * Sanity checks for the t-distribution
 */

namespace internal
{

template<typename T>
statslib_constexpr
bool
t_sanity_check(const T dof_par)
noexcept
{
    return( GCINT::is_nan(dof_par) ? \
                false :
            //
            STLIM<T>::epsilon() > dof_par ? \
                false :
            //
                true );
}

template<typename T>
statslib_constexpr
bool
t_sanity_check(const T inp_val, const T dof_par)
noexcept
{
    return (!GCINT::is_nan(inp_val)) && t_sanity_check(dof_par);
}

}
