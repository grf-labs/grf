# Test for the "bayesopt" package for auto-parameters tuning on GRF

The main function is /core/test/forest/BayesoptTest.cpp

In this test case, we first need to set the parameters for "bayesopt" package, such as number of interations, kernel name, etc. Please see line 120-138 in "BayesoptTest.cpp". And the introduction of these parameters can be found in [Bayesopt Parameters](https://rmcantin.bitbucket.io/html/usemanual.html#params)

The second step is to set up the parameters we need to tune, such as mtry, min_node_size. And set the boundary for these parameters. See line 148-153.

The third step is to write a evaluation function or loss function for our grf, please see the function testFunction_test, which calcuates the MSE of random forest for prediction.

BTW, in the CMakeList.txt file, we need to include the "bayesopt" package like "TARGET_LINK_LIBRARIES(grf bayesopt ${CMAKE_SOURCE_DIR}/bayesopt/lib/libnlopt.a)"