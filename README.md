[![Build Status](https://travis-ci.org/imbs-hl/ranger.svg?branch=master)](https://travis-ci.org/imbs-hl/ranger)
[![Coverage Status](https://coveralls.io/repos/github/imbs-hl/ranger/badge.svg?branch=master)](https://coveralls.io/github/imbs-hl/ranger?branch=master)
![CRAN Downloads month](http://cranlogs.r-pkg.org/badges/ranger?color=brightgreen)
![CRAN Downloads overall](http://cranlogs.r-pkg.org/badges/grand-total/ranger?color=brightgreen)
## ranger: A Fast Implementation of Random Forests
Marvin N. Wright, wright@imbs.uni-luebeck.de

### Introduction
Ranger is a fast implementation of random forest (Breiman 2001) or recursive partitioning, particularly suited for high dimensional data. Classification, regression, probability estimation and survival forests are supported. Classification and regression forests are implemented as in the original Random Forest (Breiman 2001), survival forests as in Random Survival Forests (Ishwaran et al. 2008). For probability estimation forests see Malley et al. (2012). 

Ranger is written in C++, but a version for R is available, too. We recommend to use the R version. It is easy to install and use and the results are readily available for further analysis. The R version is as fast as the pure C++ version.

### Installation
To install the Ranger R package from CRAN, just run

```R
install.packages("ranger”)
```

R version >= 3.1 is required. Note that, for now, R-devel and the new RTools toolchain is required for multithreading on Windows platforms (or install a binary version).

To install the C++ version of Ranger in Linux or Mac OS X you will need a compiler supporting C++11 (i.e. gcc >= 4.7 or Clang >= 3.0) and Cmake. To build start a terminal from the Ranger main directory and run the following commands

```bash
cd source
mkdir build
cd build
cmake ..
make
```

After compilation there should be an executable called "ranger" in the build directory. 

To run the C++ version in Microsoft Windows please cross compile or ask for a binary.

### Usage
For usage of the R version see ?ranger in R. Most importantly, see the Examples section. As a first example you could try 

```R  
ranger(Species ~ ., data = iris)
```

In the C++ version type 

```bash
ranger --help 
```

for a list of commands. First you need a training dataset in a file. This file should contain one header line with variable names and one line with variable values per sample. Variable names must not contain any whitespace, comma or semicolon. Values can be seperated by whitespace, comma or semicolon but can not be mixed in one file. A typical call of Ranger would be for example

```bash
ranger --verbose --file data.dat --depvarname Species --treetype 1 --ntree 1000 --nthreads 4
```

If you find any bugs, or if you experience any crashes, please report to us. If you have any questions just ask, we won't bite. 

### References
* Wright, M. N. & Ziegler, A. (2016). ranger: A Fast Implementation of Random Forests for High Dimensional Data in C++ and R. Journal of Statistical Software, in press. http://arxiv.org/abs/1508.04409.
* Schmid, M., Wright, M. N. & Ziegler, A. (2015). On the Use of Harrell's C for Node Splitting in Random Survival Forests. Technical Report. http://arxiv.org/abs/1507.03092.
* Wright, M. N., Dankowski, T. & Ziegler, A. (2016). Random forests for survival analysis using maximally selected rank statistics. Technical Report. http://arxiv.org/abs/1605.03391.
* Breiman, L. (2001). Random forests. Machine learning, 45(1), 5-32.
* Ishwaran, H., Kogalur, U. B., Blackstone, E. H., & Lauer, M. S. (2008). Random survival forests. The Annals of Applied Statistics, 841-860.
* Malley, J. D., Kruppa, J., Dasgupta, A., Malley, K. G., & Ziegler, A. (2012). Probability machines: consistent probability estimation using nonparametric learning machines. Methods Inf Med, 51(1), 74.
