#include "testfunctions.hpp"
#include "bopt_state.hpp"

/*
 * Preconditions: test_save.exe test was called before this one
 *  and generated state.dat file is accesible
 */
int main()
{
    std::cout << "IMPORTANT NOTE: Remember to execute \"test_save.exe\" before this one" << std::endl;
    std::cout << "Restoring Optimization..." << std::endl;
    
    // Second optimization (restored from first optimization state)
    bayesopt::Parameters par2;
    par2.n_iterations = 190;
    par2.random_seed = 0;
    par2.verbose_level = 1;
    par2.noise = 1e-10;

    BraninNormalized branin2(par2);

    // Restore operation and run optimization
    bayesopt::BOptState state2;
    state2.loadFromFile("state.dat", par2);
    branin2.restoreOptimization(state2);
    for(size_t i = branin2.getCurrentIter(); i < par2.n_iterations; i++){
        branin2.stepOptimization();
    }
    
    vectord result = branin2.getFinalResult();
    std::cout << "Branin2 Result: " << result << "->" 
        << branin2.evaluateSample(result) << std::endl;
        
    // Try to remove used .dat file
    if( remove( "state.dat" ) == 0 ){
        std::cout << "File \"state.dat\" successfully removed" << std::endl;
    }
    else{
        std::cout << "Error: cannot remove \"state.dat\" file" << std::endl; 
    }
    
    return 0;
}
