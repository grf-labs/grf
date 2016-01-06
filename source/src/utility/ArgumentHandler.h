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

#ifndef ARGUMENTHANDLER_H_
#define ARGUMENTHANDLER_H_

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
  MemoryMode memmode;
  bool savemem;
  std::string predict;
  std::string splitweights;
  uint nthreads;

  // All command line arguments as member: Small letters
  bool predall;
  std::vector<std::string> catvars;
  std::string file;
  ImportanceMode impmeasure;
  uint targetpartitionsize;
  uint mtry;
  std::string outprefix;
  bool probability;
  SplitRule splitrule;
  std::string statusvarname;
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

#endif /* ARGUMENTHANDLER_H_ */
