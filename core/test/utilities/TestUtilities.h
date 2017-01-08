#ifndef GRADIENTFOREST_TESTUTILITIES_H
#define GRADIENTFOREST_TESTUTILITIES_H

#include "globals.h"
#include "Observations.h"

class TestUtilities {
public:
    static Observations create_observations(std::vector<double> outcome);
    static Observations create_observations(std::vector<double> outcome,
                                            std::vector<double> treatment,
                                            std::vector<double> instrument);
};


#endif //GRADIENTFOREST_TESTUTILITIES_H
