# Author: Denis-Michael Lux
# Date: 11. January 2015
#
# This software may be modified and distributed under the terms
# of the BSD license. See the LICENSE file for details.

####################################################################
##  COMPILATION SECTION                                           ##
####################################################################
CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)

SET(UZLMATH_MAJOR 0)
SET(UZLMATH_MINOR 1)
SET(UZLMATH_PATCH 1)
SET(UZLMATH_VERSION ${UZLMATH_MAJOR}.${UZLMATH_MINOR}.${UZLMATH_PATCH})

MESSAGE(STATUS "Configuring UZLMathLib ${UZLMATH_VERSION}")

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_aux/Modules/")

PROJECT(UZLMathLib CXX)

MESSAGE(STATUS "")
MESSAGE(STATUS "CMAKE_SYSTEM_NAME          = ${CMAKE_SYSTEM_NAME}"         )
MESSAGE(STATUS "CMAKE_CXX_COMPILER_ID      = ${CMAKE_CXX_COMPILER_ID}"     )
MESSAGE(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
MESSAGE(STATUS "CMAKE_COMPILER_IS_GNUCXX   = ${CMAKE_COMPILER_IS_GNUCXX}"  )
MESSAGE(STATUS "")

MESSAGE(STATUS "*** Try to find OpenBLAS")
# Including OpenBLAS
INCLUDE(UZLMATH_FindOpenBLAS)
MESSAGE(STATUS "~> OpenBLAS_FOUND = ${OpenBLAS_FOUND}")

IF(OpenBLAS_FOUND)
    SET(UZLMATH_LIBS ${UZLMATH_LIBS} ${OpenBLAS_LIB})
    SET(UZLMATH_INCLUDE_DIRS ${UZLMATH_INCLUDE_DIRS} ${OpenBLAS_INCLUDE_DIR})

    MESSAGE(STATUS "")
    MESSAGE(STATUS "*** If the OpenBLAS library is installed in")
    MESSAGE(STATUS "*** /usr/local/lib or /usr/local/lib64")
    MESSAGE(STATUS "*** make sure the run-time linker can find it.")
    MESSAGE(STATUS "*** On Linux systems this can be done by editing /etc/ld.so.conf")
    MESSAGE(STATUS "*** or modifying the LD_LIBRARY_PATH environment variable.")
    MESSAGE(STATUS "")
ENDIF()

MESSAGE(STATUS "*** Try to find FFTW")
# Including FFTW
INCLUDE(UZLMATH_FindFFTW)
MESSAGE(STATUS "~> FFTW_FOUND = ${FFTW_FOUND}")

IF(FFTW_FOUND)
    SET(UZLMATH_LIBS ${UZLMATH_LIBS} ${FFTW_LIB})
    SET(UZLMATH_INCLUDE_DIRS ${UZLMATH_INCLUDE_DIRS} ${FFTW_INCLUDE_DIR})
ENDIF()

MESSAGE(STATUS "")
MESSAGE(STATUS "*** UZLMath wrapper library will use the following libraries:")
MESSAGE(STATUS "*** UZLMATH_LIBS          = ${UZLMATH_LIBS}"                  )
MESSAGE(STATUS "*** UZLMATH_INCLUDE_DIRS  = ${UZLMATH_INCLUDE_DIRS}"          )
MESSAGE(STATUS "")

IF(DEFINED CMAKE_CXX_COMPILER_ID AND DEFINED CMAKE_CXX_COMPILER_VERSION)
    IF(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND NOT ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 4.8.1)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
        MESSAGE(STATUS "Detected gcc 4.8.1 or later. Added '-std=c++11 -Wall -fopenmp' to compiler flags")
    ENDIF()
ENDIF()

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}/include                 )
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}/include/uzlmath_headers )
INCLUDE_DIRECTORIES( ${UZLMATH_INCLUDE_DIRS}                             )

IF(APPLE)
    IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
        #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weverything")
        MESSAGE(STATUS "Clang compiler on MacOS X detected. Added '-stdlib=libc++ -std=c++11 -O3' to compiler flags")
    ENDIF()
ENDIF()

IF(CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
ENDIF()

MESSAGE(STATUS "")
MESSAGE(STATUS "CMAKE_CXX_FLAGS           = ${CMAKE_CXX_FLAGS}"          )
MESSAGE(STATUS "CMAKE_SHARED_LINKER_FLAGS = ${CMAKE_SHARED_LINKER_FLAGS}")
MESSAGE(STATUS "CMAKE_REQUIRED_INCLUDES   = ${CMAKE_REQUIRED_INCLUDES}"  )
MESSAGE(STATUS "")

# define source files
FILE(GLOB SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.c ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp )
FILE(GLOB SRCS ${SRCS} ${CMAKE_CURRENT_SOURCE_DIR}/include/uzlmath_headers/*.h ${CMAKE_CURRENT_SOURCE_DIR}/include/uzlmath_headers/*.hpp )
FILE(GLOB SRCS ${SRCS} ${CMAKE_CURRENT_SOURCE_DIR}/include/uzlmath )

####################################################################
##  COMPILING LIBRARY                                        	  ##
####################################################################
ADD_LIBRARY(           UZLMath SHARED ${SRCS}                                                   )
TARGET_LINK_LIBRARIES( UZLMath ${UZLMATH_LIBS}                                                  )
SET_TARGET_PROPERTIES( UZLMath PROPERTIES VERSION ${UZLMATH_VERSION} SOVERSION ${UZLMATH_MAJOR} )

####################################################################
##  BENCHMARKS                                        		      ##
####################################################################
ADD_SUBDIRECTORY( benchmark )

####################################################################
##  UNIT TESTING                                        		  ##
####################################################################
ADD_SUBDIRECTORY( gtest-1.7.0   )
ADD_SUBDIRECTORY( test          )


# Creating binary and linking library
ADD_EXECUTABLE(
    soft_test
    ${PROJECT_SOURCE_DIR}/soft/soft_test.cpp
)
TARGET_LINK_LIBRARIES( soft_test UZLMath )


####################################################################
##  INSTALLATION SECTION                                          ##
####################################################################
