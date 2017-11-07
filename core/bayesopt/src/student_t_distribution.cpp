#include "student_t_distribution.hpp"


namespace bayesopt
{

  StudentTDistribution::StudentTDistribution(randEngine& eng): 
    ProbabilityDistribution(eng), d_(2)
  {
    mean_ = 0.0;  std_ = 1.0; dof_ = 2;
  }

  StudentTDistribution::~StudentTDistribution(){}


  double StudentTDistribution::negativeExpectedImprovement(double min,
							   size_t g)
  {
    const double diff = min - mean_;
    const double z = diff / std_;
  
    assert((g == 1) && "Students t EI with exponent not yet supported.");
    return -(diff * boost::math::cdf(d_,z) 
	     + (dof_*std_+z*diff)/(dof_-1) * boost::math::pdf(d_,z) ); 
  }  // negativeExpectedImprovement


  double StudentTDistribution::lowerConfidenceBound(double beta)
  {    
    return mean_ - beta*std_/sqrt(static_cast<double>(dof_));
  }  // lowerConfidenceBound

  double StudentTDistribution::negativeProbabilityOfImprovement(double min,
								double epsilon)
  {  
    return -cdf(d_,(min - mean_ + epsilon)/std_);
  }  // negativeProbabilityOfImprovement


  double StudentTDistribution::sample_query()
  { 
    double n = static_cast<double>(dof_);
    randNFloat normal(mtRandom,normalDist(mean_,std_));
    randGFloat gamma(mtRandom,gammaDist(n/2.0));
    return normal() / sqrt(2*gamma()/n);
  }  // sample_query

} //namespace bayesopt
