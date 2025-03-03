##########################################
# CMake Updated: 19th March, 2018
# 
##########################################

cmake_minimum_required (VERSION 2.6)

#For Clang
#set(CMAKE_CXX_COMPILER /sw/summitdev/llvm/20170502/clang/20170502/bin/clang)

#For XLC
#set(CMAKE_CXX_COMPILER /sw/summitdev/xl/20161123/xlC/13.1.5/bin/xlc++)
##set(CMAKE_CXX_COMPILER /sw/summitdev/xl/20170606-beta/xlC/14.1.0/bin/xlc++)

#For GCC:
#set(CMAKE_CXX_COMPILER /sw/summitdev/gcc/6.3.1-20170301/bin/g++)
#set(CMAKE_CXX_COMPILER /sw/summitdev/gcc/5.4.0new/bin/g++)

#project name
project (apply_Htarget)

#C++ 11 support -- Must have cmake > 3.0
set (CMAKE_CXX_STANDARD 11)

#For XLC
#set (CMAKE_CXX_FLAGS "-std=gnu++11")

#Add compiler flags 
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

#Adding OpenMP

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if(IBM_SUMMIT)
#Adding BLAS ESSL (Use either essl/esslsmp/esslsmpcuda)
INCLUDE_DIRECTORIES( $ENV{OLCF_ESSL_ROOT}/include )
	FIND_LIBRARY( ESSL_LIB NAMES essl PATHS "$ENV{OLCF_ESSL_ROOT}/lib64" NO_DEFAULT_PATH )
	IF ( NOT ESSL_LIB )
    	MESSAGE(FATAL_ERROR "ESSL LIB NOT FOUND")
	ENDIF()
	#Map the BLAS library path to ESSL LIB
	set (BLAS_LIB ESSL_LIB)
	MESSAGE( "Using ESSL" )
endif()

#Enable this if cmake version is less than 3.0
#IF(UNIX)
#    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x")
#ENDIF()

IF(ARM)
	#TBD - fix the search for ARMPL lib
	FIND_LIBRARY( ARMPL_LIB NAMES libarmpl armpl HINTS $ENV{ARMPL_LIBRARIES} PATHS $ENV{ARMPL_LIBRARIES} NO_DEFAULT_PATH )
	MESSAGE ("ARMPL LIB is  ${ARMPL_LIB}" ) 
	set (BLAS_LIB ${ARMPL_LIB} )
	MESSAGE ("BLASL LIB is  ${BLAS_LIB}" ) 
ENDIF()

#Enable this if running on MacOS
#IF(APPLE) 
#  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -framework IOKit -framework Cocoa -framework Foundation")
#  SET(LINK_FLAGS "${LINK_FLAGS} -framework IOKit -framework Cocoa -framework Foundation")
#  add_definitions(-DDARWIN)
#ENDIF()

find_package( BLAS REQUIRED )

#Include directories for all the include files
include_directories(${PROJECT_SOURCE_DIR}
                    ${PROJECT_SOURCE_DIR}/util
                    ${PROJECT_SOURCE_DIR}/Analysis)

#Add the two libraries as sub-directories. The CMakeLists files in these directories ensure they get compiled
add_subdirectory (util)
add_subdirectory (Analysis)


#Executable
add_executable(test_traversal ${PROJECT_SOURCE_DIR}/traverse.cpp
                              ${PROJECT_SOURCE_DIR}/apply_Htarget.cpp 
                              ${PROJECT_SOURCE_DIR}/Matrix.cpp)

add_executable(simulateFake ${PROJECT_SOURCE_DIR}/simulate_driver_generate.cpp
                            ${PROJECT_SOURCE_DIR}/apply_Htarget.cpp
                            ${PROJECT_SOURCE_DIR}/Matrix.cpp)

#Link Libraries

#test_traversal
target_link_libraries(test_traversal Util)
target_link_libraries(test_traversal ${BLAS_LIB})
target_link_libraries(test_traversal ${BLAS_LIBRARIES})
target_link_options(test_traversal PRIVATE ${BLAS_LINKER_FLAGS})

#simulateFake
target_link_libraries(simulateFake Util)
target_link_libraries(simulateFake Analysis)
target_link_libraries(simulateFake ${BLAS_LIB})
target_link_libraries(simulateFake ${BLAS_LIBRARIES})
target_link_options(simulateFake PRIVATE ${BLAS_LINKER_FLAGS})

#Link Libraries -- Additional if C++ Experimental Filesystem with GCC
#target_link_libraries(test_traversal -lstdc++fs)
