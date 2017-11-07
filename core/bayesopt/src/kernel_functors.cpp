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

#include <stdexcept>
#include "log.hpp"
#include "parser.hpp"
#include "ublas_extra.hpp"
#include "kernel_functors.hpp"

#include "kernels/kernel_atomic.hpp"
#include "kernels/kernel_const.hpp"
#include "kernels/kernel_linear.hpp"
#include "kernels/kernel_hamming.hpp"
#include "kernels/kernel_matern.hpp"
#include "kernels/kernel_polynomial.hpp"
#include "kernels/kernel_gaussian.hpp"
#include "kernels/kernel_rq.hpp"

#include "kernels/kernel_combined.hpp"
#include "kernels/kernel_sum.hpp"
#include "kernels/kernel_prod.hpp"

namespace bayesopt
{

  KernelFactory::KernelFactory()
  {
    registry["kConst"] = & create_func<ConstKernel>;
    registry["kLinear"] = & create_func<LinKernel>;
    registry["kLinearARD"] = & create_func<LinKernelARD>;

    registry["kHamming"] = & create_func<HammingKernel>;

    registry["kMaternISO1"] = & create_func<MaternIso1>;
    registry["kMaternISO3"] = & create_func<MaternIso3>;
    registry["kMaternISO5"] = & create_func<MaternIso5>;
    registry["kMaternARD1"] = & create_func<MaternARD1>;
    registry["kMaternARD3"] = & create_func<MaternARD3>;
    registry["kMaternARD5"] = & create_func<MaternARD5>;
  
    registry["kPoly1"] = & create_func<Polynomial>;
    registry["kPoly2"] = & create_func<Polynomial2>;
    registry["kPoly3"] = & create_func<Polynomial3>;
    registry["kPoly4"] = & create_func<Polynomial4>;
    registry["kPoly5"] = & create_func<Polynomial5>;
    registry["kPoly6"] = & create_func<Polynomial6>;

    registry["kSEARD"] = & create_func<SEArd>;
    registry["kSEISO"] = & create_func<SEIso>;

    registry["kRQISO"] = & create_func<RQIso>;

    registry["kSum"] = & create_func<KernelSum>;
    registry["kProd"] = & create_func<KernelProd>;
  }


  /** 
   * \brief Factory model for kernel functions
   * This function is based on the libgp library by Manuel Blum
   *      https://bitbucket.org/mblum/libgp
   * which follows the squeme of GPML by Rasmussen and Nickisch
   *     http://www.gaussianprocess.org/gpml/code/matlab/doc/
   * @param name string with the kernel structure
   * @param imput_dim number of input dimensions
   * @return kernel pointer
   */
  Kernel* KernelFactory::create(std::string name, size_t input_dim)
  {
    Kernel *kFunc;
    std::string os, os1, os2;
    utils::parseExpresion(name,os,os1,os2);

    std::map<std::string,KernelFactory::create_func_definition>::iterator it = registry.find(os);
    if (it == registry.end()) 
      {
	throw std::invalid_argument("Error while parsing kernel function: "
				    "Kernel not found " + os);
	return NULL;
      } 
    kFunc = it->second();
    if (os1.length() == 0 && os2.length() == 0) 
      {
	kFunc->init(input_dim);
      } 
    else // Combined kernel
      {
	kFunc->init(input_dim, create(os1,input_dim), create(os2,input_dim));
      }
    return kFunc;
  };


  //////////////////////////////////////////////////////////////////////

  KernelModel::KernelModel(size_t dim, Parameters parameters)
  { setKernel(parameters.kernel,dim);  }

  void KernelModel::setKernel (const vectord &thetav, 
			      const vectord &stheta,
			      std::string k_name, 
			      size_t dim)
  {
    KernelFactory mKFactory;

    mKernel.reset(mKFactory.create(k_name, dim));

    if ((thetav.size() == 1) && (stheta.size() == 1) && (mKernel->nHyperParameters() != 1))
      {
	// We assume isotropic prior, so we replicate the vectors for all dimensions
	size_t n = mKernel->nHyperParameters();

	FILE_LOG(logINFO) << "Expected " << n << " hyperparameters."
			  << " Replicating parameters and prior.";

	vectord newthetav = svectord(n,thetav(0));
	vectord newstheta = svectord(n,stheta(0));

	setKernelPrior(newthetav,newstheta);
	mKernel->setHyperParameters(newthetav);
      }
    else
      {
	setKernelPrior(thetav,stheta);
	mKernel->setHyperParameters(thetav);
      }
  }

  void KernelModel::setKernel (KernelParameters kernel, 
			      size_t dim)
  {
    size_t n = kernel.hp_mean.size();
    vectord th = kernel.hp_mean;
    vectord sth = kernel.hp_std;
    setKernel(th, sth, kernel.name, dim);
  };


  void KernelModel::computeCorrMatrix(const vecOfvec& XX, matrixd& corrMatrix, 
				     double nugget)
  {
    assert(corrMatrix.size1() == XX.size());
    assert(corrMatrix.size2() == XX.size());
    const size_t nSamples = XX.size();
  
    for (size_t ii=0; ii< nSamples; ++ii)
      {
	for (size_t jj=0; jj < ii; ++jj)
	  {
	    corrMatrix(ii,jj) = (*mKernel)(XX[ii], XX[jj]);
	    corrMatrix(jj,ii) = corrMatrix(ii,jj);
	  }
	corrMatrix(ii,ii) = (*mKernel)(XX[ii],XX[ii]) + nugget;
      }
  }

  void KernelModel::computeDerivativeCorrMatrix(const vecOfvec& XX, 
					       matrixd& corrMatrix,
					       int dth_index)
  {
    assert(corrMatrix.size1() == XX.size());
    assert(corrMatrix.size2() == XX.size());
    const size_t nSamples = XX.size();
   
    for (size_t ii=0; ii< nSamples; ++ii)
      {
	for (size_t jj=0; jj < ii; ++jj)
	  {
	    corrMatrix(ii,jj) = mKernel->gradient(XX[ii],XX[jj], 
						  dth_index);
	    corrMatrix(jj,ii) = corrMatrix(ii,jj);
	  }
	corrMatrix(ii,ii) = mKernel->gradient(XX[ii],XX[ii],dth_index);
      }
  }

  
  double KernelModel::kernelLogPrior()
  {
    double prior = 0.0;
    vectord th = mKernel->getHyperParameters();
    for(size_t i = 0; i<th.size();++i)
      {
	if (priorKernel[i].standard_deviation() > 0)
	  {
	    prior += std::log(boost::math::pdf(priorKernel[i],th(i)));
	  }
      }
    return prior;
  }
} //namespace bayesopt
