//
// Created by KuangKun on 11/3/17.
//

#ifndef GRF_BAYESOPTSRC_H
#define GRF_BAYESOPTSRC_H

#include <bayesopt/bayesopt.hpp>              // For the C++ API
#include <bayesopt/parameters.hpp>              // For the C++ API
#include <boost/numeric/ublas/assignment.hpp> // <<= op assigment
#include "commons/Data.h"

class BayesoptSrc: public bayesopt::ContinuousModel
{
private:
    Data* data;
    uint outcome_index;
public:

    BayesoptSrc(size_t dim,bayesopt::Parameters param);

    double evaluateSample( const vectord &Xi );

    //double evaluateSample_with_5000_trees( const vectord &Xi );

    bool checkReachability( const vectord &query );

    void setData(Data* data_temp, uint outcome_index_temp);

};


#endif //GRF_BAYESOPTSRC_H
