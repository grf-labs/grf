#include "DataDouble.h"

DataDouble::DataDouble() :
    data(0) {
}

DataDouble::~DataDouble() {
  if (!externalData) {
    delete[] data;
  }
}

