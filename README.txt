

ranger: A Fast Implementation of Random Forests

Written by:

  Marvin N. Wright
  Institut für Medizinische Biometrie und Statistik
  Universität zu Lübeck
  Ratzeburger Allee 160
  23562 Lübeck

  http://www.imbs-luebeck.de
  wright@imbs.uni-luebeck.de

1. Introduction

  Ranger is a fast implementation of Random Forest (Breiman 2001) or recursive partitioning, particularly suited for high dimensional data. Classification, regression, and survival forests are supported. Classification and regression forests are implemented as in the original Random Forest (Breiman 2001), survival forests as in Random Survival Forests (Ishwaran et al. 2008). 

  Ranger is written in C++, but a version for R is avaiable, too. We recommend to use the R version. It is easy to install and use and the results are readily available for further analysis. The R version is as fast as the pure C++ version.

2. Installation

  To install the Ranger R package from CRAN, just run

    install.packages("ranger”) 

  R version >= 3.1 is required. Note that, for now, no multithreading is supported in the R version on Windows platforms (the compiler in RTools is too old).

  To install the C++ version of Ranger in Linux or Mac OS X you will need a compiler supporting C++11 (i.e. gcc >= 4.7 or Clang >= 3.0) and Cmake. To build start a terminal from the Ranger main directory and run the following commands

    cd source
    mkdir build
    cd build
    cmake ..
    make

  After compilation there should be an executable called "ranger" in the build directory. 

  To run the C++ version in Microsoft Windows please cross compile or ask for a binary.

3. Usage

  For usage of the R version see ?ranger in R. Most importantly, see the Examples section. As a first example you could try 
  
    ranger(Species ~ ., data = iris)

  In the C++ version type 

    ranger --help 

  for a list of commands. First you need a training dataset in a file. This file should contain one header line with variable names and one line with variable values per sample. Variable names must not contain any whitespace, comma or semicolon. Values can be seperated by whitespace, comma or semicolon but can not be mixed in one file. A typical call of Ranger would be for example

    ranger --verbose --file data.dat --depvarname Species --treetype 1 --ntree 1000 --nthreads 4

  If you find any bugs, or if you experience any crashes, please report to us. If you have any questions just ask, we won't bite. 

References

  Breiman, L. (2001). Random forests. Machine learning, 45(1), 5-32.
  Ishwaran, H., Kogalur, U. B., Blackstone, E. H., & Lauer, M. S. (2008). Random survival forests. The Annals of Applied Statistics, 841-860.
