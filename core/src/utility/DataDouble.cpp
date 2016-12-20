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

#include "DataDouble.h"

DataDouble::DataDouble() :
    data(0) {
}

DataDouble::~DataDouble() {
  if (!externalData) {
    delete[] data;
  }
}

