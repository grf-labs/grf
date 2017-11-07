CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#Test for random number generator
ADD_EXECUTABLE(randtest ./testrand.cpp)

#Test for parsers
ADD_EXECUTABLE(parsetest ../utils/parser.cpp ./testparser.cpp)

#Test for grid
ADD_EXECUTABLE(gridtest ./testgrid.cpp)

#Test for MCMC
ADD_EXECUTABLE(mcmctest ../src/mcmc_sampler.cpp ./testmcmc.cpp)

#Test for Cholesky decomposition
ADD_EXECUTABLE(choltest ./testchol.cpp)

#Test Restore
ADD_EXECUTABLE(test_restore ./test_restore.cpp )
add_dependencies(test_restore bayesopt)
TARGET_LINK_LIBRARIES(test_restore bayesopt)

#Test Save
ADD_EXECUTABLE(test_save ./test_save.cpp )
add_dependencies(test_save bayesopt)
TARGET_LINK_LIBRARIES(test_save bayesopt)

#Test File Parser
ADD_EXECUTABLE(test_fileparser ./test_fileparser.cpp)
add_dependencies(test_fileparser bayesopt)
TARGET_LINK_LIBRARIES(test_fileparser bayesopt)

#Test Initial Samples
ADD_EXECUTABLE(test_initial_samples ./test_initial_samples.cpp)
add_dependencies(test_initial_samples bayesopt)
TARGET_LINK_LIBRARIES(test_initial_samples bayesopt)
