//
// Created by KuangKun on 11/3/17.
//

#include "BayesoptSrc.h"

#include "commons/utility.h"
#include "forest/ForestPredictor.h"
#include "forest/ForestPredictors.h"
#include "forest/ForestTrainer.h"
#include "forest/ForestTrainers.h"

void init_trainer_test(ForestTrainer& trainer,
                       bool honesty,
                       uint ci_group_size,
                       uint number_of_tree,
                       const double *query,
                       uint variable_dimensions) {
//uint mtry = 3;
    uint mtry = query[0]*variable_dimensions;
    uint num_trees = number_of_tree;
    uint seed = 42;
    uint num_threads = 4;
//uint min_node_size = 1;
    uint min_node_size = query[1];
    std::set<size_t> no_split_variables;
    std::string split_select_weights_file = "";
    bool sample_with_replacement = true;
    std::string sample_weights_file = "";
//double sample_fraction = ci_group_size > 1 ? 0.35 : 0.7;
    double sample_fraction = query[2];

    trainer.init(mtry, num_trees, seed, num_threads,
                 min_node_size, no_split_variables, split_select_weights_file,
                 sample_with_replacement, sample_weights_file, sample_fraction,
                 honesty, ci_group_size);
}

double testFunction_test(Data* data, uint outcome_index, unsigned int n, const double *query, uint number_of_trees,
                         double *gradient, /* NULL if not needed */
                         void *func_data)
{
    /*
    char buffer[100];
    getcwd(buffer, 100);
    std::cout << "The   current   directory   is: " << buffer << std::endl;
    //*/

    //Data* data = load_data("/Users/kuangkun/Documents/Stanford/Generalized Random Forests/Code/grftuning_bak/core/test/forest/dataset/AirFoil.txt");
    //uint outcome_index = data->get_num_cols()-1;
    uint variable_dimensions = data->get_num_cols()-1;


    double alpha = 0.10;

    //uint mtry = query[0];
    //uint min_node_size = query[1];

    ForestTrainer trainer = ForestTrainers::regression_trainer(data, outcome_index, alpha);
    init_trainer_test(trainer, false, 2, number_of_trees, query, variable_dimensions);


    Forest forest = trainer.train(data);
    ForestPredictor predictor = ForestPredictors::regression_predictor(4, 2);
    std::vector<Prediction> predictions = predictor.predict_oob(forest, data);
    double MSE = 0;
    for (int i = 0; i < predictions.size(); ++i){
        double y_real = data->get(i,outcome_index);
        double y_pred = predictions[i].get_predictions().at(0);
        //std::cout << "*" << i << y_real-y_pred << std::endl;
        MSE += (y_real-y_pred)*(y_real-y_pred);
    }
    //delete  data;

    MSE = MSE / predictions.size();
    return MSE;
}

BayesoptSrc::BayesoptSrc(size_t dim, bayesopt::Parameters param):
    ContinuousModel(dim,param) {};

double BayesoptSrc::evaluateSample( const vectord &Xi )
{
    double x[100];
    for (size_t i = 0; i < Xi.size(); ++i)
    {
        x[i] = Xi(i);
    }
    //return testFunction(Xi.size(),x,NULL,NULL);
    return testFunction_test(data, outcome_index, Xi.size(),x,20,NULL,NULL);
};

/*
double BayesoptSrc::evaluateSample_with_5000_trees( const vectord &Xi )
{
    double x[100];
    for (size_t i = 0; i < Xi.size(); ++i)
    {
        x[i] = Xi(i);
    }
    //return testFunction(Xi.size(),x,NULL,NULL);
    return testFunction_test(Xi.size(),x,5000,NULL,NULL);
};
//*/


bool BayesoptSrc::checkReachability( const vectord &query )
{ return true; };


void BayesoptSrc::setData(Data* data_temp, uint outcome_index_temp){

    data = data_temp;
    outcome_index = outcome_index_temp;
}


