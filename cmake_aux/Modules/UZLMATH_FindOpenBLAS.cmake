SET(Open_BLAS_INCLUDE_SEARCH_PATHS
  	/home/617286/Schreibtisch/OpenBLAS/include
	/usr/include
	/usr/include/openblas-base
	/usr/local/include
	/usr/local/include/openblas-base
	/usr/local/Cellar/openblas/0.2.14/include
  	/opt/OpenBLAS/include
  	$ENV{OpenBLAS_HOME}
  	$ENV{OpenBLAS_HOME}/include
)

SET(Open_BLAS_LIB_SEARCH_PATHS
	/home/617286/Schreibtisch/OpenBLAS/lib
	/lib/
	/lib/openblas-base
	/lib64/
	/usr/lib
	/usr/lib/openblas-base
	/usr/lib64
	/usr/local/lib
	/usr/local/lib64
	/usr/local/Cellar/openblas/0.2.14/lib
	/opt/OpenBLAS/lib
	$ENV{OpenBLAS}
	$ENV{OpenBLAS}/lib
	$ENV{OpenBLAS_HOME}
	$ENV{OpenBLAS_HOME}/lib
)

FIND_PATH(OpenBLAS_INCLUDE_DIR NAMES cblas.h PATHS ${Open_BLAS_INCLUDE_SEARCH_PATHS})
FIND_LIBRARY(OpenBLAS_LIB NAMES openblas PATHS ${Open_BLAS_LIB_SEARCH_PATHS})

SET(OpenBLAS_FOUND TRUE)

# Check include files
IF(NOT OpenBLAS_INCLUDE_DIR)
    SET(OpenBLAS_FOUND FALSE)
    MESSAGE(STATUS "~> Could not find OpenBLAS include. Switch OpenBLAS_FOUND to FALSE")
    IF(APPLE)
    	MESSAGE(STATUS "~> To install OpenBLAS you can use 'Homebrew'")
    	MESSAGE(STATUS "~> More information about Homebrew can be found at")
    	MESSAGE(STATUS "~>     'http://brew.sh'")
    	MESSAGE(STATUS "~> To install OpenBLAS via Homebrew run the following")
    	MESSAGE(STATUS "~> command on your Terminal")
    	MESSAGE(STATUS "~>     'brew install homebrew/science/openblas'")
    ELSE()
    	MESSAGE(STATUS "~> To install OpenBLAS get the source code by")
    	MESSAGE(STATUS "~> cloning the repository from github.")
    	MESSAGE(STATUS "~> ")
    	MESSAGE(STATUS "~> Open your Terminal and run")
    	MESSAGE(STATUS "~>     'git clone git@github.com:xianyi/OpenBLAS'")
    	MESSAGE(STATUS "~> After cloning the repository follow the instructions on")
    	MESSAGE(STATUS "~>     https://github.com/xianyi/OpenBLAS/wiki/Installation-Guide")
    ENDIF()
ENDIF()

# Check libraries
IF(NOT OpenBLAS_LIB)
    SET(OpenBLAS_FOUND FALSE)
    MESSAGE(STATUS "~> Could not find OpenBLAS lib. Switch OpenBLAS_FOUND to FALSE")
ENDIF()

IF(OpenBLAS_FOUND)
	IF (NOT OpenBLAS_FIND_QUIETLY)
    	MESSAGE(STATUS "~> Found OpenBLAS libraries: ${OpenBLAS_LIB}"        )
		MESSAGE(STATUS "~> Found OpenBLAS include:   ${OpenBLAS_INCLUDE_DIR}")
	ENDIF (NOT OpenBLAS_FIND_QUIETLY)
ELSE(OpenBLAS_FOUND)
	IF (OpenBLAS_FIND_REQUIRED)
    	MESSAGE(FATAL_ERROR "~> Could not find OpenBLAS")
	ENDIF (OpenBLAS_FIND_REQUIRED)
ENDIF(OpenBLAS_FOUND)

MARK_AS_ADVANCED(
    OpenBLAS_INCLUDE_DIR
    OpenBLAS_LIB
    OpenBLAS
)