#include "log.hpp"
#include "parser.hpp"
#include "ublas_extra.hpp"
#include "mean_functors.hpp"
#include "mean_atomic.hpp"
#include "mean_combined.hpp"

namespace bayesopt
{
  
  MeanFactory::MeanFactory()
  {
    registry["mZero"] = & create_func<ZeroFunction>;
    registry["mOne"] = & create_func<OneFunction>;
    registry["mConst"] = & create_func<ConstantFunction>;
    registry["mLinear"] = & create_func<LinearFunction>;
    registry["mSum"] = & create_func<SumFunction>;
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
  ParametricFunction* MeanFactory::create(std::string name, size_t input_dim)
  {
    ParametricFunction *mFunc;
    std::string os, os1, os2;
    utils::parseExpresion(name,os,os1,os2);

    std::map<std::string,MeanFactory::create_func_definition>::iterator it = registry.find(os);
    if (it == registry.end()) 
      {
	FILE_LOG(logERROR) << "Error: Fatal error while parsing mean function: "
			   << os << " not found" << std::endl;
	return NULL;
      } 
    mFunc = registry.find(os)->second();
    if (os1.length() == 0 && os2.length() == 0) 
      {
	mFunc->init(input_dim);
      } 
    else 
      {
	mFunc->init(input_dim, create(os1,input_dim), 
		    create(os2,input_dim));
      }
    return mFunc;

  };

  //////////////////////////////////////////////////////////////////////

  MeanModel::MeanModel(size_t dim, Parameters parameters)
  {
    setMean(parameters.mean,dim);
  }

  void MeanModel::setMean (const vectord &muv,
				     const vectord &smu,
				     std::string m_name,
				     size_t dim)
  {
    MeanFactory mPFactory;

    mMean.reset(mPFactory.create(m_name,dim));

    if ("mZero" == m_name) 
      {
	mMu = zvectord(1);
	mS_Mu = svectord(1,1e-10);
      }
    else if("mOne" == m_name) 
      {
	mMu = svectord(1,1.0);
	mS_Mu = svectord(1,1e-10);
      }
    else
      {
	mMu = muv; mS_Mu = smu;
      }

    mMean->setParameters(mMu);
  }

  void MeanModel::setMean (MeanParameters mean, size_t dim)
  {
    size_t n_mu = mean.coef_mean.size();
    vectord vmu = mean.coef_mean;
    vectord smu = mean.coef_std;
    setMean(vmu, smu, mean.name, dim);
  };

}//namespace bayesopt



