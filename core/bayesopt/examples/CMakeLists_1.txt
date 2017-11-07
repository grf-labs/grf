CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#Optimization with system calls
ADD_EXECUTABLE(branin_system_calls ./branin_system_calls.cpp)
add_dependencies(branin_system_calls bayesopt)
TARGET_LINK_LIBRARIES(branin_system_calls bayesopt)

#Optimization with xml and system calls
ADD_EXECUTABLE(branin_xml ./branin_xml.cpp)
add_dependencies(branin_xml bayesopt)
TARGET_LINK_LIBRARIES(branin_xml bayesopt)

#Continuous test
ADD_EXECUTABLE(bo_cont ./bo_cont.cpp)
add_dependencies(bo_cont bayesopt)
TARGET_LINK_LIBRARIES(bo_cont bayesopt)

#Discrete test
ADD_EXECUTABLE(bo_disc ./bo_disc.cpp)
add_dependencies(bo_disc bayesopt)
TARGET_LINK_LIBRARIES(bo_disc bayesopt)

#1D test
ADD_EXECUTABLE(bo_oned ./bo_oned.cpp)
add_dependencies(bo_oned bayesopt)
TARGET_LINK_LIBRARIES(bo_oned bayesopt)

#Branin
ADD_EXECUTABLE(bo_branin ./bo_branin.cpp )
add_dependencies(bo_branin bayesopt)
TARGET_LINK_LIBRARIES(bo_branin bayesopt)

ADD_EXECUTABLE(bo_branin_timed ./bo_branin_timed.cpp )
add_dependencies(bo_branin_timed bayesopt)
TARGET_LINK_LIBRARIES(bo_branin_timed bayesopt)

ADD_EXECUTABLE(bo_branin_mcmc ./bo_branin_mcmc.cpp )
add_dependencies(bo_branin_mcmc bayesopt)
TARGET_LINK_LIBRARIES(bo_branin_mcmc bayesopt)

#Hartmann
ADD_EXECUTABLE(bo_hartmann ./bo_hartmann.cpp )
add_dependencies(bo_hartmann bayesopt)
TARGET_LINK_LIBRARIES(bo_hartmann bayesopt)

#Camelback
ADD_EXECUTABLE(bo_camelback ./bo_camelback.cpp )
add_dependencies(bo_camelback bayesopt)
TARGET_LINK_LIBRARIES(bo_camelback bayesopt)

ADD_EXECUTABLE(bo_compare ./bo_compare.cpp )
add_dependencies(bo_compare bayesopt)
TARGET_LINK_LIBRARIES(bo_compare bayesopt)

enable_language(C)  #For OpenGL and bo_display demo

#Display test
find_package(GLUT)
find_package(OpenGL)
if(OPENGL_FOUND AND GLUT_FOUND)

  SET(DISPLAY_SRC     
    ../utils/displaygp.cpp
    )

  if(GLUT_LIBRARY MATCHES freeglut_static.lib)
      add_definitions(-DFREEGLUT_STATIC)
  endif()


  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/matplotpp
		      ${GLUT_INCLUDE_DIRS} 
		      ${GLUT_INCLUDE_DIR} 
		      ${OpenGL_INCLUDE_DIRS})			  

  link_directories(${GLUT_LIBRARY_DIRS}
		   ${OpenGL_LIBRARY_DIRS})

  add_definitions(${GLUT_DEFINITIONS}
		  ${OpenGL_DEFINITIONS})

  ADD_LIBRARY(matplotpp STATIC 
    ../matplotpp/matplotpp.cc 
    ../matplotpp/gl2ps.c
    )

  TARGET_LINK_LIBRARIES(matplotpp
    ${GLUT_LIBRARY} ${OPENGL_LIBRARY})

  ADD_EXECUTABLE(bo_display 
    ./bo_display.cpp 
    ${DISPLAY_SRC})
  add_dependencies(bo_display bayesopt matplotpp)
  TARGET_LINK_LIBRARIES(bo_display bayesopt matplotpp
    ${GLUT_LIBRARY} ${OPENGL_LIBRARY})

  ADD_EXECUTABLE(bo_branin_display
    ./bo_branin_display.cpp 
    ${DISPLAY_SRC})
  add_dependencies(bo_branin_display bayesopt matplotpp)
  TARGET_LINK_LIBRARIES(bo_branin_display bayesopt matplotpp
    ${GLUT_LIBRARY} ${OPENGL_LIBRARY})


endif()
