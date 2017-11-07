# -------------------------------------------------------------------------
#    This file is part of BayesOpt, an efficient C++ library for 
#    Bayesian optimization.
#
#    Copyright (C) 2011-2015 Ruben Martinez-Cantin <rmcantin@unizar.es>
# 
#    BayesOpt is free software: you can redistribute it and/or modify it 
#    under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    BayesOpt is distributed in the hope that it will be useful, but 
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with BayesOpt.  If not, see <http://www.gnu.org/licenses/>.
# ------------------------------------------------------------------------

PROJECT(BayesOpt CXX)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} 
		      ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(APPLE)
  # For macports or similar systems
  set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /opt/local/lib /opt/local/Library)
  set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} /opt/local/include
                                               /opt/local/Library)
  set(CMAKE_PROGRAM_PATH ${CMAKE_PROGRAM_PATH} /opt/local/bin/ /opt/local/Library)
endif(APPLE)


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING 
      "Choose the type of build, options are: Debug Release 
      RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

option(BAYESOPT_BUILD_EXAMPLES "Build examples and demos?" OFF)
option(BAYESOPT_BUILD_TESTS "Build tests?" OFF)
option(BAYESOPT_PYTHON_INTERFACE "Build Python interface?" OFF)
option(BAYESOPT_MATLAB_COMPATIBLE "Build library compatible with Matlab?" OFF)
option(BAYESOPT_BUILD_SOBOL "Build support for Sobol sequences?" ON)
option(BAYESOPT_BUILD_SHARED "Build BayesOpt as a shared library?" ON)

find_package( Boost REQUIRED )
if(Boost_FOUND)
   include_directories(${Boost_INCLUDE_DIRS})
else(Boost_FOUND)
   find_library(Boost boost PATHS /opt/local/lib)
   include_directories(${Boost_LIBRARY_PATH})
endif()

FIND_LIBRARY(NLOPT nlopt)
IF(NLOPT MATCHES NLOPT-NOTFOUND)
  SET(NLOPT_BUILD ON CACHE BOOL "Build included version of NLOPT?")
ELSE(NLOPT MATCHES NLOPT-NOTFOUND)
  SET(NLOPT_BUILD OFF CACHE BOOL "Build included version of NLOPT?")
ENDIF(NLOPT MATCHES NLOPT-NOTFOUND)

INCLUDE(UseDoxygen)

# Sobol sequences are hardcoded tables, so it might take a lot of time
# to compile.
IF(BAYESOPT_BUILD_SOBOL)
  ADD_DEFINITIONS(-DUSE_SOBOL)	
  SET(SOBOL_SRC
    ./sobol/sobol.cpp
    ./sobol/sobol_i4.cpp
    ./sobol/sobol_i8.cpp
    )
ELSE(BAYESOPT_BUILD_SOBOL)
  SET(SOBOL_SRC )
ENDIF(BAYESOPT_BUILD_SOBOL)


SET( BAYESOPT_SRCS
  ./src/bayesoptcont.cpp
  ./src/bayesoptdisc.cpp
  ./src/bayesoptbase.cpp
  ./src/bopt_state.cpp
  ./src/posteriormodel.cpp
  ./src/posterior_fixed.cpp
  ./src/posterior_empirical.cpp
  ./src/posterior_mcmc.cpp
  ./src/mcmc_sampler.cpp
  ./src/inneroptimization.cpp
  ./src/dataset.cpp
  ./src/nonparametricprocess.cpp
  ./src/kernelregressor.cpp
  ./src/conditionalbayesprocess.cpp
  ./src/gaussian_process.cpp
  ./src/gaussian_process_ml.cpp
  ./src/gaussian_process_normal.cpp
  ./src/gaussian_process_hierarchical.cpp
  ./src/student_t_process_jef.cpp
  ./src/student_t_process_nig.cpp
  ./src/parameters.cpp
  ./src/kernel_functors.cpp
  ./src/criteria_functors.cpp
  ./src/criteria_hedge.cpp
  ./src/mean_functors.cpp
  ./src/gauss_distribution.cpp
  ./src/student_t_distribution.cpp
  )

SET(UTILS_SRC
  ./utils/fileparser.cpp
  ./utils/param_loader.cpp
  ./utils/parser.cpp
  ./utils/ublas_extra.cpp
  ${SOBOL_SRC}
  )

SET(WRAPPPERS_SRC 
  ./src/wrappers/bayesoptwpr.cpp
  )


INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/include 
		     ${PROJECT_SOURCE_DIR}/utils
		     ${PROJECT_SOURCE_DIR}/sobol )


LINK_DIRECTORIES( ${PROJECT_BINARY_DIR}/lib )

IF(NLOPT_BUILD)
  ADD_SUBDIRECTORY(nlopt2)
  include_directories(${PROJECT_SOURCE_DIR}/nlopt2/api)
  SET(EXT_LIBS nlopt)
ELSE(NLOPT_BUILD)
  SET(EXT_LIBS ${NLOPT})
ENDIF(NLOPT_BUILD)


IF(BAYESOPT_BUILD_SHARED)
  ADD_LIBRARY(bayesopt SHARED ${BAYESOPT_SRCS}   
    ${WRAPPPERS_SRC} ${UTILS_SRC} )
  IF(WIN32)
    ADD_DEFINITIONS(-DBAYESOPT_DLL)	
    # In new versions of CMAKE they use a different system and the
    # symbol is not defined
    ADD_DEFINITIONS(-Dbayesopt_EXPORT )
  ENDIF()
ELSE()
  ADD_LIBRARY(bayesopt STATIC ${BAYESOPT_SRCS}   
    ${WRAPPPERS_SRC} ${UTILS_SRC} )
ENDIF()


IF(NLOPT_BUILD)
  add_dependencies(bayesopt nlopt)
ENDIF(NLOPT_BUILD)


IF((BAYESOPT_BUILD_SHARED OR BAYESOPT_MATLAB_COMPATIBLE) AND NOT WIN32)
  SET_TARGET_PROPERTIES(bayesopt PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF()
  
TARGET_LINK_LIBRARIES(bayesopt ${EXT_LIBS})

IF(BAYESOPT_BUILD_TESTS)
  ADD_SUBDIRECTORY(tests)
endif(BAYESOPT_BUILD_TESTS)


#IF(BAYESOPT_BUILD_EXAMPLES)
  #ADD_EXECUTABLE(bo_branin ./examples/bo_branin.cpp )
  #add_dependencies(bo_branin bayesopt)
  #TARGET_LINK_LIBRARIES(bo_branin bayesopt)
  #Continuous test
  #ADD_EXECUTABLE(bo_cont ./examples/bo_cont.cpp)
  #add_dependencies(bo_cont bayesopt)
  #TARGET_LINK_LIBRARIES(bo_cont bayesopt)
  
  #ADD_SUBDIRECTORY(examples)
#endif(BAYESOPT_BUILD_EXAMPLES)

# INSTALL(FILES 
#   ./include/bayesopt.hpp
#   ./include/bayesoptbase.hpp
#   ./include/boparameters.hpp
#   ./include/boparameters.h
#   ./wrappers/bayesopt.h 
#   DESTINATION include
# )

INSTALL(
  TARGETS bayesopt
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

INSTALL(
  DIRECTORY include/bayesopt DESTINATION include
)


IF(BAYESOPT_PYTHON_INTERFACE)
  INCLUDE(PythonMagic)
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
#  find_package(python)
#  find_python_module(numpy REQUIRED)
#  string(REGEX REPLACE "__init.py.*" "core/include/${module}" 
#	    		PY_MOD_INCLUDE "PY_NUMPY")
#  include_directories(${PYTHON_INCLUDE_DIRS} ${PY_NUMPY_INCLUDE_DIRS} )

  ADD_LIBRARY(esopt MODULE ${BAYESOPT_SRCS}   
    ${WRAPPPERS_SRC} ${UTILS_SRC} ./python/bayesopt.cpp)

  IF(NLOPT_BUILD)
    add_dependencies(esopt nlopt)
  ENDIF(NLOPT_BUILD)

  TARGET_LINK_LIBRARIES(esopt ${EXT_LIBS} ${PYTHON_LIBRARIES} )

  IF(NOT WIN32)
    # Kind of a hack but it works
    SET_TARGET_PROPERTIES(esopt PROPERTIES PREFIX "bay" SUFFIX ".so")
  ELSE()
    SET_TARGET_PROPERTIES(esopt PROPERTIES PREFIX "libbay" SUFFIX ".dll")
  ENDIF()

  #Find where to install Python libs
  execute_process ( COMMAND 
         python -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
         OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)

  INSTALL(
    TARGETS esopt
    LIBRARY DESTINATION ${PYTHON_SITE_PACKAGES}
    ARCHIVE DESTINATION ${PYTHON_SITE_PACKAGES}
  )
  INSTALL(
    FILES ./python/bayesoptmodule.py
    DESTINATION ${PYTHON_SITE_PACKAGES}
  )

ENDIF(BAYESOPT_PYTHON_INTERFACE)
