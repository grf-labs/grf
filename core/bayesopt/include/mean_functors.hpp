/** \file mean_functors.hpp \brief Mean (parametric) functions. */
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

#ifndef  _MEAN_FUNCTORS_HPP_
#define  _MEAN_FUNCTORS_HPP_

#include <map>
#include <boost/scoped_ptr.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "bayesopt/parameters.hpp"
#include "specialtypes.hpp"

namespace bayesopt
{
  
  /**\addtogroup ParametricFunctions
   *  \brief Set of parametric models for surrogate modelling
   */
  //@{

  /** \brief Interface for mean functors */
  class ParametricFunction
  {
  public:
    virtual int init(size_t input_dim) {return 0;};
    virtual int init(size_t input_dim, ParametricFunction* left, 
		     ParametricFunction* right) {return 0;};

    virtual void setParameters(const vectord& params) = 0;
    virtual vectord getParameters() = 0;
    virtual size_t nParameters() = 0;

    virtual double getMean(const vectord& x) = 0;
    virtual vectord operator()(const vecOfvec& x)
    {
      vectord result(x.size());
      std::vector<vectord>::const_iterator x_it = x.begin();
      std::vector<vectord>::const_iterator x_end = x.end();
      vectord::iterator res_it = result.begin();
      while(x_it != x_end)
	{
	  *res_it++ = getMean(*x_it++);
	}
      return result;
    };


    virtual size_t nFeatures() = 0;
    virtual vectord getFeatures(const vectord& x) = 0;
    virtual matrixd getAllFeatures(const vecOfvec& x)
    {
      size_t nf = nFeatures();
      matrixd result(nf,x.size());

      for(size_t ii=0; ii< x.size(); ++ii)
	{
	  column(result,ii) = getFeatures(x[ii]);
	}

      return result;
    };
  
    virtual ~ParametricFunction(){};

  protected:
    size_t n_inputs;
  };


  template <typename MeanType> ParametricFunction * create_func()
  {
    return new MeanType();
  }


  /** 
   * \brief Factory model for parametric functions
   * This factory is based on the libgp library by Manuel Blum
   *      https://bitbucket.org/mblum/libgp
   * which follows the squeme of GPML by Rasmussen and Nickisch
   *     http://www.gaussianprocess.org/gpml/code/matlab/doc/
   */
  class MeanFactory
  {
  public:
    MeanFactory ();
    virtual ~MeanFactory () {};
  
    ParametricFunction* create(std::string name, size_t input_dim);
    
  private:
    typedef ParametricFunction* (*create_func_definition)();
    std::map<std::string , MeanFactory::create_func_definition> registry;
  };

  class MeanModel
  {
  public:
    MeanModel(size_t dim, Parameters parameters);
    virtual ~MeanModel() {};

    ParametricFunction* getMeanFunc();
   
    void setParameters(const vectord &theta);
    vectord getParameters();
    size_t nParameters();

    vectord getFeatures(const vectord& x);  
    void getFeatures(const vectord& x, vectord& kx);  
    size_t nFeatures();

    void setPoints(const vecOfvec &x);
    void addNewPoint(const vectord &x);

    vectord muTimesFeat();
    double muTimesFeat(const vectord& x);

    /** 
     * \brief Select the parametric part of the surrogate process.
     * 
     * @param muv mean function parameters
     * @param smu std function parameters
     * @param m_name mean function name
     * @param dim number of input dimensions
     * @return error_code
     */
    void setMean (const vectord &muv, const vectord &smu, 
		 std::string m_name, size_t dim);

    /** Wrapper of setMean for the C++ structure */
    void setMean (MeanParameters mean, size_t dim);

    matrixd mFeatM;           ///< Value of the mean features at the input points

  private:
    vectord mMu;                 ///< Mean of the parameters of the mean function
    vectord mS_Mu;    ///< Variance of the params of the mean function W=mS_Mu*I

    boost::scoped_ptr<ParametricFunction> mMean;    ///< Pointer to mean function   
  };



  inline ParametricFunction* MeanModel::getMeanFunc()
  { return mMean.get(); }

  inline void MeanModel::setParameters(const vectord &theta)
  { mMean->setParameters(theta); };
    
  inline vectord MeanModel::getParameters(){return mMean->getParameters();};
  inline size_t MeanModel::nParameters(){return mMean->nParameters();};

  inline vectord MeanModel::getFeatures(const vectord& x) 
  { return mMean->getFeatures(x); }

  inline void MeanModel::getFeatures(const vectord& x, vectord& kx)
  { kx = mMean->getFeatures(x); }  

  inline size_t MeanModel::nFeatures()
  { return mMean->nFeatures(); }  

  inline void MeanModel::setPoints(const vecOfvec &x)
  { mFeatM = mMean->getAllFeatures(x); };

  inline void MeanModel::addNewPoint(const vectord &x)
  { 
    using boost::numeric::ublas::column;
    
    mFeatM.resize(mFeatM.size1(),mFeatM.size2()+1);  
    column(mFeatM,mFeatM.size2()-1) = mMean->getFeatures(x);
  }

  inline vectord MeanModel::muTimesFeat()
  {  return boost::numeric::ublas::prod(mMu,mFeatM); }
    
  inline double MeanModel::muTimesFeat(const vectord& x)
  { return boost::numeric::ublas::inner_prod(mMu,mMean->getFeatures(x));}



  //@}

} //namespace bayesopt


#endif
