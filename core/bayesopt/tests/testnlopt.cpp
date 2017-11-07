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


/**
 * Unit testing for the NLOPT wrappers
 * 
 */

#include <iostream>
#include <vector>
#include "FastDelegate.h"

class FooInterface
{
public:
  FooInterface(){i=0.0;}
  virtual ~FooInterface(){};
  virtual double call(const double& step)
  { 
    i+=step; 
    //    std::cout << "Foo Interface: " << i << std::endl;
    return i;
  }
protected:
  double i;
};

class FooDerivate: public FooInterface
{
public:
  FooDerivate(){j=5.0;}
  virtual ~FooDerivate(){};
  double call(const double& step)
  {
    j=step/j;
    //    std::cout << "Foo Derivate: " << j << std::endl;
    return j;
  }
private:
  double j;  
};


class Optimizer
{
public:
  Optimizer(){iter = 10e6;}
  virtual ~Optimizer(){};
  void setCallback(FooInterface* foo)
  {foo_ = static_cast<void*>(foo);}

  static double step(void* ptr, const double& query)
  {
    FooInterface* p = static_cast<FooInterface*>(ptr);
    return p->call(query);
  }
  void run()
  {
    std::vector<double> res;
    for(size_t i=10;i<iter;++i)
      {
	res.push_back(step(foo_,static_cast<double>(i)));
      }
  }
private:
  size_t iter;
  void* foo_;
};

typedef fastdelegate::FastDelegate1<const double&, double> myfunc;

class OptimizerDelegate
{
public:
  OptimizerDelegate(){iter = 10e6;}
  virtual ~OptimizerDelegate(){};
  void setCallback(FooInterface* foo)
  {f_.bind(foo, &FooInterface::call);}
  void run()
  {
    std::vector<double> res;
    for(size_t i=10;i<iter;++i)
      {
	res.push_back(f_(static_cast<double>(i)));
      }
  }

private:
  size_t iter;
  myfunc f_;
};


//////////////////////////////////////////////////////////////////////

class OptimizerBase
{
public:
  OptimizerBase(){iter = 10e6;}
  virtual ~OptimizerBase(){};
  virtual double call(const double& step) = 0;
  void run()
  {
    std::vector<double> res;
    for(size_t i=10;i<iter;++i)
      {
	res.push_back(call(static_cast<double>(i)));
      }
  }
private:
  size_t iter;
};



class BarInterface: public OptimizerBase
{
public:
  BarInterface(){i=0.0;}
  virtual ~BarInterface(){};
  virtual double call(const double& step)
  { 
    i+=step; 
    //    std::cout << "Bar Interface: " << i << std::endl;
    return i;
  }
protected:
  double i;
};

class BarDerivate: public BarInterface
{
public:
  BarDerivate(){j=5.0;}
  virtual ~BarDerivate(){};
  double call(const double& step)
  {
    j=step/j;
    // std::cout << "Bar Derivate: " << j << std::endl;
    return j;
  }
private:
  double j;  
};


int main( int argc, const char* argv[] )
{
  FooInterface* ptr1 = new FooDerivate();
  BarInterface* ptr2 = new BarDerivate();
  FooInterface* ptr3 = new FooDerivate();

  Optimizer op1;
  op1.setCallback(ptr1);

  Optimizer op2;
  op2.setCallback(ptr3);

  clock_t start = clock();
  op1.run();
  double dif1 = (double)(clock()-start) / (double)CLOCKS_PER_SEC;

  start = clock();
  ptr2->run();
  double dif2 = (double)(clock()-start) / (double)CLOCKS_PER_SEC;

  start = clock();
  op2.run();
  double dif3 = (double)(clock()-start) / (double)CLOCKS_PER_SEC;

  std::cout << dif1 << "|||" << dif2 << "|||" << dif3 << std::endl;
  return 0;
}

