/** \file kernel_functors.hpp \brief Kernel (covariance) functions */
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

#ifndef  _KERNEL_FUNCTORS_HPP_
#define  _KERNEL_FUNCTORS_HPP_

#include <map>
#include <boost/scoped_ptr.hpp>
#include <boost/math/distributions/normal.hpp> 
#include "bayesopt/parameters.hpp"
#include "specialtypes.hpp"

namespace bayesopt
{
  
  /**\addtogroup KernelFunctions
   * \brief Set of kernel or covariance functions for the nonparametric
   * processes.
   */
  //@{

  /** \brief Interface for kernel functors */
  class Kernel
  {
  public:
    virtual ~Kernel(){};
    virtual void init(size_t input_dim) {};
    virtual void init(size_t input_dim, Kernel* left, Kernel* right) {};

    virtual void setHyperParameters(const vectord &theta) = 0;
    virtual vectord getHyperParameters() = 0;
    virtual size_t nHyperParameters() = 0;

    virtual double operator()( const vectord &x1, const vectord &x2 ) = 0;
    virtual double gradient( const vectord &x1, const vectord &x2,
			     size_t component ) = 0;

  protected:
    size_t n_inputs;
  };



  template <typename KernelType> Kernel * create_func()
  {
    return new KernelType();
  }

  /** 
   * \brief Factory model for kernel functions
   * This factory is based on the libgp library by Manuel Blum
   *      https://bitbucket.org/mblum/libgp
   * which follows the squeme of GPML by Rasmussen and Nickisch
   *     http://www.gaussianprocess.org/gpml/code/matlab/doc/
   */
  class KernelFactory
  {
  public:
    KernelFactory ();
    virtual ~KernelFactory () {};
  
    Kernel* create(std::string name, size_t input_dim);
    
  private:
    typedef Kernel* (*create_func_definition)();
    std::map<std::string , KernelFactory::create_func_definition> registry;
  };


  class KernelModel
  {
  public:
    KernelModel(size_t dim, Parameters parameters);
    virtual ~KernelModel() {};

    Kernel* getKernel();
    
    void setHyperParameters(const vectord &theta);
    vectord getHyperParameters();
    size_t nHyperParameters();


    /** 
     * \brief Select kernel (covariance function) for the surrogate process.
     * @param thetav kernel parameters (mean)
     * @param stheta kernel parameters (std)
     * @param k_name kernel name
     */
    void setKernel (const vectord &thetav, const vectord &stheta, 
		   std::string k_name, size_t dim);

    /** Wrapper of setKernel for C++ kernel structure */
    void setKernel (KernelParameters kernel, size_t dim);

    void computeCorrMatrix(const vecOfvec& XX, matrixd& corrMatrix, double nugget);
    void computeDerivativeCorrMatrix(const vecOfvec& XX, matrixd& corrMatrix, 
				    int dth_index);
    vectord computeCrossCorrelation(const vecOfvec& XX, const vectord &query);
    void computeCrossCorrelation(const vecOfvec& XX, const vectord &query,
				 vectord& knx);
    double computeSelfCorrelation(const vectord& query);
    double kernelLogPrior();

  private:
    /** Set prior (Gaussian) for kernel hyperparameters */
    void setKernelPrior (const vectord &theta, const vectord &s_theta);

    boost::scoped_ptr<Kernel> mKernel;            ///< Pointer to kernel function
    std::vector<boost::math::normal> priorKernel; ///< Prior of kernel parameters
  };

  inline Kernel* KernelModel::getKernel()
  { return mKernel.get();  }

  inline void KernelModel::setHyperParameters(const vectord &theta)
  { mKernel->setHyperParameters(theta); };
    
  inline vectord KernelModel::getHyperParameters()
  {return mKernel->getHyperParameters();};
  
  inline size_t KernelModel::nHyperParameters()
  {return mKernel->nHyperParameters();};

  inline vectord KernelModel::computeCrossCorrelation(const vecOfvec& XX, 
						      const vectord &query)
  {
    vectord knx(XX.size());
    computeCrossCorrelation(XX,query,knx);
    return knx;
  }

  inline void KernelModel::computeCrossCorrelation(const vecOfvec& XX, 
						   const vectord &query,
						   vectord& knx)
  {
    std::vector<vectord>::const_iterator x_it  = XX.begin();
    vectord::iterator k_it = knx.begin();
    while(x_it != XX.end())
      {	*k_it++ = (*mKernel)(*x_it++, query); }
  }


  inline double KernelModel::computeSelfCorrelation(const vectord& query)
  { return (*mKernel)(query,query); }

  inline void KernelModel::setKernelPrior (const vectord &theta, 
					   const vectord &s_theta)
  {
    for (size_t i = 0; i<theta.size(); ++i)
      {
	boost::math::normal n(theta(i),s_theta(i));
	priorKernel.push_back(n);
      }
  };


  //@}

} //namespace bayesopt


#endif
