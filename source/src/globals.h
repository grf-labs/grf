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

#ifndef GLOBALS_H_
#define GLOBALS_H_

#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
    TypeName(const TypeName&);             \
    void operator=(const TypeName&)

typedef unsigned int uint;

// Tree types, probability is not selected by ID
enum TreeType {
  TREE_CLASSIFICATION = 1,
  TREE_REGRESSION = 3,
  TREE_SURVIVAL = 5,
  TREE_PROBABILITY = 9
};

// Memory modes
enum MemoryMode {
  MEM_DOUBLE = 0,
  MEM_FLOAT = 1,
  MEM_CHAR = 2
};
const uint MAX_MEM_MODE = 2;

// Mask and Offset to store 2 bit values in bytes
static const int mask[4] = {192,48,12,3};
static const int offset[4] = {6,4,2,0};

// Variable importance
enum ImportanceMode {
  IMP_NONE = 0,
  IMP_GINI = 1,
  IMP_PERM_BREIMAN = 2,
  IMP_PERM_LIAW = 4,
  IMP_PERM_RAW = 3
};
const uint MAX_IMP_MODE = 4;

// Default values
const uint DEFAULT_NUM_TREE = 500;
const uint DEFAULT_NUM_THREADS = 0;
const ImportanceMode DEFAULT_IMPORTANCE_MODE = IMP_NONE;

const uint DEFAULT_MIN_NODE_SIZE_CLASSIFICATION = 1;
const uint DEFAULT_MIN_NODE_SIZE_REGRESSION = 5;
const uint DEFAULT_MIN_NODE_SIZE_SURVIVAL = 3;
const uint DEFAULT_MIN_NODE_SIZE_PROBABILITY = 10;

// Interval to print progress in seconds
const double STATUS_INTERVAL = 30.0;

// Threshold for q value split method switch
const double Q_THRESHOLD = 0.02;

#endif /* GLOBALS_H_ */
