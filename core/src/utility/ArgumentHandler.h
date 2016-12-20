#ifndef GRADIENTFOREST_ARGUMENTHANDLER_H_
#define GRADIENTFOREST_ARGUMENTHANDLER_H_

#include <getopt.h>
#include <string>
#include <vector>

#include "globals.h"

/*
 * Encapsulate getopt.
 * To add an option:
 *    Add to short_options
 *    Add to long_options
 *    Add member variable
 *    Add default value to constructor
 *    Add case in processArguments() function, use try-catch
 *    Add to checkArguments() function?
 *    Add to help function
 *    Add in R version?
 * Access via public members
 */
class ArgumentHandler {
public:
  ArgumentHandler(int argc, char **argv);
  virtual ~ArgumentHandler();

  // Get arguments and catch conversion exceptions
  int processArguments();

  // Check required arguments, ranges, files, ..
  void checkArguments();

  // All command line arguments as member: Capital letters
  std::vector<std::string> alwayssplitvars;
  std::string caseweights;
  std::string depvarname;
  double fraction;
  bool savemem;
  std::string predict;
  std::string splitweights;
  uint nthreads;

  // All command line arguments as member: Small letters
  std::string file;
  uint targetpartitionsize;
  uint mtry;
  std::vector<double>* quantiles;
  std::string statusvarname;
  std::string instrumentvarname;
  uint ntree;
  bool replace;
  bool verbose;
  bool write;
  TreeType treetype;
  uint seed;

private:
  // Display messages
  void displayHelp();
  void displayVersion();

  int argc;
  char** argv;

  DISALLOW_COPY_AND_ASSIGN(ArgumentHandler);
};

#endif /* GRADIENTFOREST_ARGUMENTHANDLER_H_ */
