BayesOpt: A Bayesian optimization library
=========================================

BayesOpt is an efficient implementation of the Bayesian optimization
methodology for nonlinear optimization, experimental design and
hyperparameter tunning.

Bayesian optimization uses a distribution over functions to build a
surrogate model of the unknown function for we are looking the
optimum, and then apply some active learning strategy to select the
query points that provides most potential interest or
improvement. Thus, it is a **sample efficient** method for nonlinear
optimization, design of experiments and simulations or bandits-like
problems. Currently, it is being used in many scientific and
industrial applications. In the literature it is also called
Sequential Kriging Optimization (SKO), Sequential Model-Based
Optimization (SMBO) or Efficient Global Optimization (EGO).

BayesOpt is licensed under the AGPL and it is free to use. However,
if you use BayesOpt in a work that leads to a scientific
publication, we would appreciate it if you would kindly **cite BayesOpt**
in your manuscript.

> Ruben Martinez-Cantin, **BayesOpt: A Bayesian Optimization
> Library for Nonlinear Optimization, Experimental Design and
> Bandits**. Journal of Machine Learning Research, 15(Nov):3735--3739, 2014.

The paper can be found at http://jmlr.org/papers/v15/martinezcantin14a.html

Commercial applications may also acquire a **commercial license**. Please
contact <rmcantin@unizar.es> for details.


Getting and installing BayesOpt
-------------------------------

The library can be download from any of this sources:

- Download: <https://bitbucket.org/rmcantin/bayesopt>
- Mirror: <https://github.com/rmcantin/bayesopt>
- Mirror: <http://mloss.org/software/view/453/>

You can also get the *cutting-edge* version from the repositories:

    >> hg clone https://bitbucket.org/rmcantin/bayesopt

or the git mirror:

    >> git clone https://github.com/rmcantin/bayesopt


The online documentation can be found at:
<http://rmcantin.bitbucket.io/html/> where it includes a [install
guide](http://rmcantin.bitbucket.io/html/install.html).


Questions and issues
--------------------
- The best place to ask questions and discuss about BayesOpt is the
[bayesopt-discussion mailing
list](https://groups.google.com/forum/#!forum/bayesopt-discussion). 
- Please file bug reports or suggestions at:
https://bitbucket.org/rmcantin/bayesopt/issues or https://github.com/rmcantin/bayesopt/issues
- Alternatively, you may directly contact Ruben Martinez-Cantin <rmcantin@unizar.es>.


----------------------------------------------------------------------

Copyright (C) 2011-2015 Ruben Martinez-Cantin <rmcantin@unizar.es>

BayesOpt is free software: you can redistribute it and/or modify it
under the terms of the GNU Affero General Public License as published by 
the Free Software Foundation, either version 3 of the License, or 
(at your option) any later version.

BayesOpt is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with BayesOpt. If not, see <http://www.gnu.org/licenses/>.

----------------------------------------------------------------------
