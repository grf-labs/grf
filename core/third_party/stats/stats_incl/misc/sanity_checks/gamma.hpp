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
 * Sanity checks for the Gamma distribution
 */

namespace internal
{

template<typename T>
statslib_constexpr
bool
gamma_sanity_check(const T shape_par, const T scale_par)
noexcept
{
    return( GCINT::any_nan(shape_par,scale_par) ? \
                false :
            //
            shape_par < T(0) ? \
                false :
            //
            STLIM<T>::epsilon() > scale_par ? \
                false :
            //
                true );
}

template<typename T>
statslib_constexpr
bool
gamma_sanity_check(const T inp_val, const T shape_par, const T scale_par)
noexcept
{
    return (!GCINT::is_nan(inp_val)) && gamma_sanity_check(shape_par,scale_par);
}

}
