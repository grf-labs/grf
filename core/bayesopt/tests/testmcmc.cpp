/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2015 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   BayesOpt is free software: you can redistribute it and/or modify it 
   under the terms of the GNU Affero General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   BayesOpt is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Affero General Public License for more details.

   You should have received a copy of the GNU Affero General Public License
   along with BayesOpt.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------
*/

#define _USE_MATH_DEFINES
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include "specialtypes.hpp"
#include "ublas_cholesky.hpp"
#include "posterior_mcmc.hpp"

namespace bnu = boost::numeric::ublas;
 
int determinant_sign(const bnu::permutation_matrix<std::size_t>& pm)
{
  int pm_sign=1;
  std::size_t size = pm.size();
  for (std::size_t i = 0; i < size; ++i)
    if (i != pm(i))
      pm_sign *= -1.0; // swap_rows would swap a pair of rows here, so we change sign
  return pm_sign;
}
 
double determinant(bnu::matrix<double>& m ) {
  bnu::permutation_matrix<std::size_t> pm(m.size1());
  double det = 1.0;
  if( bnu::lu_factorize(m,pm) ) {
    det = 0.0;
  } else {
    for(int i = 0; i < m.size1(); i++) 
      det *= m(i,i); // multiply by elements on diagonal
    det = det * determinant_sign( pm );
  }
  return det;
}


double gauss(const vectord& x, const vectord& mu, const matrixd& sigma)
{
  const double tpi = boost::math::constants::two_pi<double>();
  double n = static_cast<double>(x.size());
  const vectord vd = x-mu;
  matrixd invS = sigma;
  bayesopt::utils::inverse_cholesky(sigma,invS);
  matrixd sig = sigma;

  return pow(tpi,n/2)*pow(determinant(sig),0.5)*exp(-0.5*inner_prod(vd,prod(invS,vd)));
}

class Posterior: public bayesopt::RBOptimizable
{
  double evaluate(const vectord& x)
  {
    vectord mu1(2), mu2(2), mu3(2);
    matrixd s1(2,2), s2(2,2), s3(2,2);

    mu1 <<= 0,0;
    mu2 <<= 1,1;
    mu3 <<= -4,2;
    
    s1 <<= 1,0, 
      0,1;

    s2 <<= 4,0,
      0,0.6;

    s3 <<= 4,0,
      0,0.6;
  
    return gauss(x,mu1,s1) + gauss(x,mu2,s2) + gauss(x,mu3,s3);
  }
};
  
int main()
{
  randEngine reng;
  Posterior post;
  bayesopt::MCMCSampler sampler(&post,2,reng);
  vectord x = zvectord(2);
  sampler.run(x);
  sampler.printParticles();

  return 0;
}
