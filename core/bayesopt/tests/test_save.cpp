#include "testfunctions.hpp"
#include "bopt_state.hpp"

int main()
{
    // First optimization
    bayesopt::Parameters par1;
    par1.n_iterations = 190;
    par1.random_seed = 0;
    par1.verbose_level = 1;
    par1.noise = 1e-10;

    BraninNormalized branin1(par1);
    vectord result(2);
    
    // Initialize first optimization that will stop at half run
    size_t stopAt = par1.n_iterations/2;
    branin1.initializeOptimization();
    for(size_t i = branin1.getCurrentIter(); i < stopAt; i++){
        branin1.stepOptimization();
    }
    // Save state
    bayesopt::BOptState state;
    branin1.saveOptimization(state);
    state.saveToFile("state.dat");
    std::cout << "STATE ITERS: " << state.mCurrentIter << std::endl;
    
    result = branin1.getFinalResult();
    std::cout << "Branin1 Result: " << result << "->" 
        << branin1.evaluateSample(result) << std::endl;
}
