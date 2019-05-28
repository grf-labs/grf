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
 * log-if-true function
 */

#ifndef _statslib_log_if_HPP
#define _statslib_log_if_HPP

namespace internal
{

template<typename T>
statslib_constexpr
T
log_if(const T x, const bool log_form)
noexcept
{
    return log_form ? stmath::log(x) : x;
}

template<typename T>
statslib_constexpr
T
log_zero_if(const bool log_form)
noexcept
{
    return log_form ? - STLIM<T>::infinity() : T(0);
}

template<typename T>
statslib_constexpr
T
log_one_if(const bool log_form)
noexcept
{
    return log_form ? T(0) : T(1);
}

}

#endif
