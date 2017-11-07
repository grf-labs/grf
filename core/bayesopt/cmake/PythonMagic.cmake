# Some "magic" to find python in different platforms. Most of the code got from
# the openrave simulator.
# Still needs some extra magic to work seamlessly in Windows.

include(FindPythonInterp)
## check python
find_package(PythonLibs 2) # using PYTHON_INCLUDE_PATH instead of PYTHON_INCLUDE_DIR
if( NOT PYTHON_EXECUTABLE )
  # look specifically for 2.7
  find_program(PYTHON_EXECUTABLE NAMES python2.7 python27 python python.exe PATHS [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.7\\InstallPath])
endif()

if( PYTHON_EXECUTABLE )
  # architecture independent
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib(0))"
    OUTPUT_VARIABLE _python_sitepackage OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _python_failed0)
  # architexture dependent
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib(1))"
    OUTPUT_VARIABLE _python_distpackage OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _python_failed1)
#  execute_process(
#    COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; from os.path import relpath; print(relpath(get_python_lib(1,prefix='${CMAKE_INSTALL_PREFIX}'),'${CMAKE_INSTALL_PREFIX}'))"
#    OUTPUT_VARIABLE OPENRAVE_PYTHON_INSTALL_DIR OUTPUT_STRIP_TRAILING_WHITESPACE
#    RESULT_VARIABLE _python_failed2)

  if( ${_python_failed0} EQUAL 0 AND ${_python_failed1} EQUAL 0 ) #AND ${_python_failed2} EQUAL 0 )
    if( NOT IS_DIRECTORY "${_python_distpackage}/numpy/core/include" )
      if( NOT IS_DIRECTORY "${_python_sitepackage}/numpy/core/include" )
	FIND_PACKAGE(Numpy REQUIRED)
	if(Numpy MATCHES Numpy-NOTFOUND)
          message(STATUS "cannot find python numpy dir!")
          set(PYTHON_EXECUTABLE)
	else()
	  set(PYTHON_INCLUDE_PATH ${PYTHON_INCLUDE_PATH} ${PYTHON_NUMPY_INCLUDE_DIR})	   
	endif() 
      else()
        set(PYTHON_INCLUDE_PATH ${PYTHON_INCLUDE_PATH} ${_python_sitepackage_arch}/numpy/core/include)
      endif()
    else()
      set(PYTHON_INCLUDE_PATH ${PYTHON_INCLUDE_PATH} ${_python_distpackage}/numpy/core/include)
    endif()

    # get the major.minor python version
    execute_process(
      COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print('%d.%d'%sys.version_info[0:2])"
      OUTPUT_VARIABLE _python_version OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE _python_failed)
    if( ${_python_failed} EQUAL 0 )
      string(REGEX REPLACE "[\r\n]" "" PYTHON_MAJORMINOR_VERSION "${_python_version}")
    else()
      message(STATUS "failed to get python version")
    endif()
  else()
    message(STATUS "failed to get python site-package directories via get_pytho_lib")
    set(PYTHON_EXECUTABLE)
  endif()
endif()
