
#include <cmath> // std::abs
#include "testfunctions.hpp"
#include "bopt_state.hpp"

int main()  
{
    std::string filename = "test_initial_samples.dat";
    float precision = 0.0000001f;
    int saved_samples = 0; // this parameters can be modified to check other cases
    
    /*
     * One optimization generate only the initial samples.
     * Then the evaluated values are discarded and saved into file
     */
    bayesopt::Parameters params;
    params.verbose_level = 0;
    params.n_init_samples = 5;
    params.n_iterations = 10;

    BraninNormalized branin1(params);
    branin1.initializeOptimization();
    bayesopt::BOptState state1;
    branin1.saveOptimization(state1);
    
    // Save samples for later comparison
    vectord first_mY(params.n_init_samples);
    for(size_t i=0; i<params.n_init_samples; i++)
      {
        first_mY[i] = state1.mY[i]; 
      }
    state1.mY.resize(saved_samples); // remove all evaluation of samples
    state1.saveToFile(filename);

    /*
     * Second optimization restores non-evaluated samples
     * During restoreOptimization, the mX should not change
     */
    params.load_save_flag = 1;
    params.load_filename = filename;

    BraninNormalized branin2(params);
    vectord _result(2);
    branin2.optimize(_result);
    
    bayesopt::BOptState state2;
    branin2.saveOptimization(state2);
    state2.saveToFile(filename);
    
    // Try to remove used .dat file
    if( remove( filename.c_str() ) == 0 )
      {
        std::cout << "File \"" << filename << "\" successfully removed" << std::endl;
      }
    else
      {
        std::cout << "Error: cannot remove \"" << filename << "\" file" << std::endl; 
      }
    
    int returnValue = 0;
    // The length of state1.mY must be zero (as it was cleaned):
    if(state1.mY.size() != saved_samples){
        returnValue = -1;
        std::cout << "ERROR: length of state1.mY expected " << state1.mY.size()
            << " to be equals to " << saved_samples << std::endl;
    }
    
    // The evaluated initial samples at mX and mY (as bo_branin evaluations are deterministic)
    // are expected to contain the same values:
    for(size_t i=0; i<params.n_init_samples && returnValue != -1; i++){
        for(size_t j=0; j<state1.mX[i].size() && returnValue != -1; j++){
            if(std::abs(state1.mX[i][j] - state2.mX[i][j]) > precision){
                std::cout << "ERROR: mX["<< i <<"] was " << state2.mX[i] << 
                    " expected to be equals to " << state1.mX[i] << std::endl;
                returnValue = -1;
            }
            if(std::abs(first_mY[i] - state2.mY[i]) > precision){
                std::cout << "ERROR: mY["<< i <<"] was " << state2.mY[i] <<
                    " expected to be equals to " << first_mY[i] << std::endl;
                returnValue = -1;
            }
        }
    }

    if(returnValue == 0){
        std::cout << "Tests completed without errors" << std::endl;
    }
    return returnValue;
}
