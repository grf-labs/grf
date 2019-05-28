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
 * Sample from a multinomial distribution
 */

template<typename mT, typename eT>
statslib_inline
mT
rmultinom(const mT& prob)
{
    const ullint_t n_prob = mat_ops::n_elem(prob);

    ullint_t n_j = n_prob;

    //

    mT ret(n_prob,1);
    const mT prob_csum = mat_ops::cumsum(prob);

    eT p_j = prob(0,0);
    ret(0,0) = rbinom(n_j,p_j);

    //

    ullint_t ret_sum = ret(0,0);
    
    for (ullint_t j = 1U; j < n_prob; j++)
    {
        p_j = prob(j,0) / (eT(1) - prob_csum(j-1,0));
        n_j = n_prob - ret_sum;
        
        ret(j,0) = rbinom(n_j,p_j);
        
        ret_sum += ret(j,0);
    }

    //

    return ret;
}
