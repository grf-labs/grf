/*-------------------------------------------------------------------------------
This file is part of Ranger.
    
Ranger is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Ranger is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Ranger. If not, see <http://www.gnu.org/licenses/>.

Written by: 

Marvin N. Wright
Institut für Medizinische Biometrie und Statistik
Universität zu Lübeck
Ratzeburger Allee 160
23562 Lübeck 

http://www.imbs-luebeck.de
wright@imbs.uni-luebeck.de
#-------------------------------------------------------------------------------*/

#include <fstream>
#include <iostream>
#include <stdexcept>

#include "ArgumentHandler.h"
#include "version.h"
#include "utility.h"

ArgumentHandler::ArgumentHandler(int argc, char **argv) :
    caseweights(""), depvarname(""), fraction(1), memmode(MEM_DOUBLE), savemem(false), predict(""), splitweights(""), nthreads(
        DEFAULT_NUM_THREADS), predall(false), file(""), impmeasure(DEFAULT_IMPORTANCE_MODE), targetpartitionsize(0), mtry(
        0), outprefix("ranger_out"), probability(false), splitrule(DEFAULT_SPLITRULE), statusvarname(""), ntree(
        DEFAULT_NUM_TREE), replace(true), verbose(false), write(false), treetype(TREE_CLASSIFICATION), seed(0) {
  this->argc = argc;
  this->argv = argv;
}

ArgumentHandler::~ArgumentHandler() {
}

int ArgumentHandler::processArguments() {

  // short options
  char const *short_options = "A:C:D:F:M:NP:S:U:Zac:f:hil::m:o:pr:s:t:uvwy:z:";

  // long options: longname, no/optional/required argument?, flag(not used!), shortname
    const struct option long_options[] = {

      { "alwayssplitvars",      required_argument,  0, 'A'},
      { "caseweights",          required_argument,  0, 'C'},
      { "depvarname",           required_argument,  0, 'D'},
      { "fraction",             required_argument,  0, 'F'},
      { "memmode",              required_argument,  0, 'M'},
      { "savemem",              no_argument,        0, 'N'},
      { "predict",              required_argument,  0, 'P'},
      { "splitweights",         required_argument,  0, 'S'},
      { "nthreads",             required_argument,  0, 'U'},
      { "version",              no_argument,        0, 'Z'},

      { "predall",              no_argument,  0, 'a'},
      { "catvars",              required_argument,  0, 'c'},
      { "file",                 required_argument,  0, 'f'},
      { "help",                 no_argument,        0, 'h'},
      { "impmeasure",           required_argument,  0, 'i'},
      { "targetpartitionsize",  required_argument,  0, 'l'},
      { "mtry",                 required_argument,  0, 'm'},
      { "outprefix",            required_argument,  0, 'o'},
      { "probability",          no_argument,        0, 'p'},
      { "splitrule",            required_argument,  0, 'r'},
      { "statusvarname",        required_argument,  0, 's'},
      { "ntree",                required_argument,  0, 't'},
      { "noreplace",            no_argument,        0, 'u'},
      { "verbose",              no_argument,        0, 'v'},
      { "write",                no_argument,        0, 'w'},
      { "treetype",             required_argument,  0, 'y'},
      { "seed",                 required_argument,  0, 'z'},

      { 0, 0, 0, 0}
    };

  while (1) {
    int option_index = 0;
    int c = getopt_long(argc, argv, short_options, long_options, &option_index);

    // stop if no options left
    if (c == -1) {
      break;
    }

    switch (c) {

    // upper case options
    case 'A':
      splitString(alwayssplitvars, optarg, ',');
      break;

    case 'C':
      caseweights = optarg;
      break;

    case 'D':
      depvarname = optarg;
      break;

    case 'F':
      try {
        fraction = std::stod(optarg);
        if (fraction > 1 || fraction <= 0) {
          throw std::runtime_error("");
        }
      } catch (...) {
        throw std::runtime_error(
            "Illegal argument for option 'fraction'. Please give a value in (0,1]. See '--help' for details.");
      }
      break;

    case 'M':
      try {
        memmode = (MemoryMode) std::stoi(optarg);
        if (memmode > MAX_MEM_MODE) {
          throw std::runtime_error("");
        }
      } catch (...) {
        throw std::runtime_error(
            "Illegal argument for option 'memmode'. Please give a positive integer. See '--help' for details.");
      }
      break;

    case 'N':
      savemem = true;
      break;

    case 'P':
      predict = optarg;
      break;

    case 'S':
      splitweights = optarg;
      break;

    case 'U':
      try {
        int temp = std::stoi(optarg);
        if (temp < 1) {
          throw std::runtime_error("");
        } else {
          nthreads = temp;
        }
      } catch (...) {
        throw std::runtime_error(
            "Illegal argument for option 'nthreads'. Please give a positive integer. See '--help' for details.");
      }
      break;

    case 'Z':
      displayVersion();
      return -1;
      break;

    // lower case options
    case 'a':
      predall = true;
      break;

    case 'c':
      splitString(catvars, optarg, ',');
      break;

    case 'f':
      file = optarg;
      break;

    case 'h':
      displayHelp();
      return -1;
      break;

    case 'i':
      try {
        impmeasure = (ImportanceMode) std::stoi(optarg);
        if (impmeasure > MAX_IMP_MODE) {
          throw std::runtime_error("");
        }
      } catch (...) {
        throw std::runtime_error(
            "Illegal argument for option 'impmeasure'. Please give a positive integer. See '--help' for details.");
      }
      break;

    case 'l':
      try {
        int temp = std::stoi(optarg);
        if (temp < 1) {
          throw std::runtime_error("");
        } else {
          targetpartitionsize = temp;
        }
      } catch (...) {
        throw std::runtime_error(
            "Illegal argument for option 'targetpartitionsize'. Please give a positive integer. See '--help' for details.");
      }
      break;

    case 'm':
      try {
        int temp = std::stoi(optarg);
        if (temp < 1) {
          throw std::runtime_error("");
        } else {
          mtry = temp;
        }
      } catch (...) {
        throw std::runtime_error(
            "Illegal argument for option 'mtry'. Please give a positive integer. See '--help' for details.");
      }
      break;

    case 'o':
      outprefix = optarg;
      break;

    case 'p':
      probability = true;
      break;

    case 'r':
      try {
        switch (std::stoi(optarg)) {
        case 1:
          splitrule = LOGRANK;
          break;
        case 2:
          splitrule = AUC;
          break;
        case 3:
          splitrule = AUC_IGNORE_TIES;
          break;
        default:
          throw std::runtime_error("");
          break;
        }
      } catch (...) {
        throw std::runtime_error("Illegal splitrule selected. See '--help' for details.");
      }
      break;

    case 's':
      statusvarname = optarg;
      break;

    case 't':
      try {
        int temp = std::stoi(optarg);
        if (temp < 1) {
          throw std::runtime_error("");
        } else {
          ntree = temp;
        }
      } catch (...) {
        throw std::runtime_error(
            "Illegal argument for option 'ntree'. Please give a positive integer. See '--help' for details.");
      }
      break;

    case 'u':
      replace = false;
      break;

    case 'v':
      verbose = true;
      break;

    case 'w':
      write = true;
      break;

    case 'y':
      try {
        switch (std::stoi(optarg)) {
        case 1:
          treetype = TREE_CLASSIFICATION;
          break;
        case 3:
          treetype = TREE_REGRESSION;
          break;
        case 5:
          treetype = TREE_SURVIVAL;
          break;
        default:
          throw std::runtime_error("");
          break;
        }
      } catch (...) {
        throw std::runtime_error(
            "Illegal argument for option 'treetype'. Please give a positive integer. See '--help' for details.");
      }
      break;

    case 'z':
      try {
        int temp = std::stoi(optarg);
        if (temp < 0) {
          throw std::runtime_error("");
        } else {
          seed = temp;
        }
      } catch (...) {
        throw std::runtime_error(
            "Illegal argument for option 'seed'. Please give a positive integer. See '--help' for details.");
      }
      break;

    default:
      break;

    }
  }

  // print all other parameters
  while (optind < argc) {
    std::cout << "Other parameter, not processed: " << argv[optind++] << std::endl;
  }

  return 0;
}

void ArgumentHandler::checkArguments() {

  // required arguments
  if (file.empty()) {
    throw std::runtime_error("Please specify an input filename with '--file'. See '--help' for details.");
  }
  if (predict.empty() && depvarname.empty()) {
    throw std::runtime_error("Please specify a dependent variable name with '--depvarname'. See '--help' for details.");
  }

  if (treetype == TREE_SURVIVAL && statusvarname.empty()) {
    throw std::runtime_error("Please specify a status variable name with '--statusvarname'. See '--help' for details.");
  }
  if (treetype != TREE_SURVIVAL && !statusvarname.empty()) {
    throw std::runtime_error("Option '--statusvarname' only applicable for survival forest. See '--help' for details.");
  }

  if (treetype == TREE_SURVIVAL && impmeasure == IMP_GINI) {
    throw std::runtime_error(
        "Node impurity variable importance not supported for survival forests. See '--help' for details.");
  }

  if (treetype != TREE_CLASSIFICATION && probability) {
    throw std::runtime_error("Probability estimation is only applicable to classification forests.");
  }

  // Get treetype for prediction
  if (!predict.empty()) {
    std::ifstream infile;
    infile.open(predict, std::ios::binary);
    if (!infile.good()) {
      throw std::runtime_error("Could not read from input file: " + predict + ".");
    }

    // Do not read dependent_varID, num_variables, num_trees and is_ordered_variable
    infile.seekg(2 * sizeof(size_t));
    size_t length;
    infile.read((char*) &length, sizeof(length));
    infile.seekg(4 * sizeof(size_t) + length * sizeof(bool));

    // Get treetype
    infile.read((char*) &treetype, sizeof(treetype));
    infile.close();
  }

  // Option predall only for classification and regression
  if (predall && treetype != TREE_CLASSIFICATION && treetype != TREE_REGRESSION) {
    throw std::runtime_error("Option '--predall' only available for classification and regression.");
  }

  if (predict.empty() && predall) {
    throw std::runtime_error("Option '--predall' only available in prediction mode.");
  }

  if (!alwayssplitvars.empty() && !splitweights.empty()) {
    throw std::runtime_error("Please use only one option of splitweights and alwayssplitvars.");
  }

  // Check splitrule
  if ((splitrule == AUC || splitrule == AUC_IGNORE_TIES) && treetype != TREE_SURVIVAL) {
    throw std::runtime_error("Illegal splitrule selected. See '--help' for details.");
  }

}

void ArgumentHandler::displayHelp() {
  std::cout << "Usage: " << std::endl;
  std::cout << "    " << argv[0] << " [options]" << std::endl;
  std::cout << std::endl;

  std::cout << "Options:" << std::endl;
  std::cout << "    " << "--help                        Print this help." << std::endl;
  std::cout << "    " << "--version                     Print version and citation information." << std::endl;
  std::cout << "    " << "--verbose                     Turn on verbose mode." << std::endl;
  std::cout << "    " << "--file FILE                   Filename of input data. Only numerical values are supported." << std::endl;
  std::cout << "    " << "--treetype TYPE               Set tree type to:" << std::endl;
  std::cout << "    " << "                              TYPE = 1: Classification." << std::endl;
  std::cout << "    " << "                              TYPE = 3: Regression." << std::endl;
  std::cout << "    " << "                              TYPE = 5: Survival." << std::endl;
  std::cout << "    " << "                              (Default: 1)" << std::endl;
  std::cout << "    " << "--probability                 Grow a Classification forest with probability estimation for the classes." << std::endl;
  std::cout << "    " << "                              Use in combination with --treetype 1." << std::endl;
  std::cout << "    " << "--depvarname NAME             Name of dependent variable. For survival trees this is the time variable." << std::endl;
  std::cout << "    " << "--statusvarname NAME          Name of status variable, only applicable for survival trees." << std::endl;
  std::cout << "    " << "                              Coding is 1 for event and 0 for censored." << std::endl;
  std::cout << "    " << "--ntree N                     Set number of trees to N." << std::endl;
  std::cout << "    " << "                              (Default: 500)" << std::endl;
  std::cout << "    " << "--mtry N                      Number of variables to possibly split at in each node." << std::endl;
  std::cout << "    " << "                              (Default: sqrt(p) for Classification and Survival, p/3 for Regression. " << std::endl;
  std::cout << "    " << "                               p = number of independent variables)" << std::endl;
  std::cout << "    " << "--targetpartitionsize N       Set minimal node size to N." << std::endl;
  std::cout << "    " << "                              For Classification and Regression growing is stopped if a node reaches a size smaller than N." << std::endl;
  std::cout << "    " << "                              For Survival growing is stopped if one child would reach a size smaller than N." << std::endl;
  std::cout << "    " << "                              This means nodes with size smaller N can occur for Classification and Regression." << std::endl;
  std::cout << "    " << "                              (Default: 1 for Classification, 5 for Regression, and 3 for Survival)" << std::endl;
  std::cout << "    " << "--catvars V1,V2,..            Comma separated list of names of (unordered) categorical variables. " << std::endl;
  std::cout << "    " << "                              Categorical variables must contain only positive integer values." << std::endl;
  std::cout << "    " << "--write                       Save forest to file <outprefix>.forest." << std::endl;
  std::cout << "    " << "--predict FILE                Load forest from FILE and predict with new data." << std::endl;
  std::cout << "    " << "--predall                     Return a matrix with individual predictions for each tree instead of aggregated " << std::endl;
  std::cout << "    " << "                              predictions for all trees (classification and regression only)." << std::endl;
  std::cout << "    " << "--impmeasure TYPE             Set importance mode to:" << std::endl;
  std::cout << "    " << "                              TYPE = 0: none." << std::endl;
  std::cout << "    " << "                              TYPE = 1: Node impurity: Gini for Classification, variance for Regression." << std::endl;
  std::cout << "    " << "                              TYPE = 2: Permutation importance, scaled by standard errors." << std::endl;
  std::cout << "    " << "                              TYPE = 3: Permutation importance, no scaling." << std::endl;
  std::cout << "    " << "                              (Default: 0)" << std::endl;
  std::cout << "    " << "--noreplace                   Sample without replacement." << std::endl;
  std::cout << "    " << "--fraction X                  Fraction of observations to sample. Default is 1 for sampling with replacement " << std::endl;
  std::cout << "    " << "                              and 0.632 for sampling without replacement." << std::endl;
  std::cout << "    " << "--splitrule RULE              Splitting rule:" << std::endl;
  std::cout << "    " << "                              RULE = 1: Gini for Classification, variance for Regression, logrank for Survival." << std::endl;
  std::cout << "    " << "                              RULE = 2: AUC for Survival, not available for Classification and Regression." << std::endl;
  std::cout << "    " << "                              (Default: 1)" << std::endl;
  std::cout << "    " << "--caseweights FILE            Filename of case weights file." << std::endl;
  std::cout << "    " << "--splitweights FILE           Filename of split select weights file." << std::endl;
  std::cout << "    " << "--alwayssplitvars V1,V2,..    Comma separated list of variable names to be always considered for splitting." << std::endl;
  std::cout << "    " << "--nthreads N                  Set number of parallel threads to N." << std::endl;
  std::cout << "    " << "                              (Default: Number of CPUs available)" << std::endl;
  std::cout << "    " << "--seed SEED                   Set random seed to SEED." << std::endl;
  std::cout << "    " << "                              (Default: No seed)" << std::endl;
  std::cout << "    " << "--outprefix PREFIX            Prefix for output files." << std::endl;
  std::cout << "    " << "--memmode MODE                Set memory mode to:" << std::endl;
  std::cout << "    " << "                              MODE = 0: double." << std::endl;
  std::cout << "    " << "                              MODE = 1: float." << std::endl;
  std::cout << "    " << "                              MODE = 2: char." << std::endl;
  std::cout << "    " << "                              (Default: 0)" << std::endl;
  std::cout << "    " << "--savemem                     Use memory saving (but slower) splitting mode." << std::endl;
  std::cout << std::endl;

  std::cout << "See README file for details and examples." << std::endl;
}

// TODO: Change citation info
void ArgumentHandler::displayVersion() {
  std::cout << "Ranger version: " << RANGER_VERSION << std::endl;
  std::cout << std::endl;
  std::cout << "Please cite Ranger: " << std::endl;
  std::cout << "Marvin N. Wright and .. (2014). Ranger. Journal." << std::endl;
  std::cout << std::endl;
  std::cout << "BibTeX:" << std::endl;
  std::cout << "@Article{," << std::endl;
  std::cout << "    title = {Ranger}" << std::endl;
  std::cout << "    author = {Marvin N. Wright and ..}," << std::endl;
  std::cout << "    journal = {Journal}," << std::endl;
  std::cout << "    year = {2014}," << std::endl;
  std::cout << "}" << std::endl;
}
