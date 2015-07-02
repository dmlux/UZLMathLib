# Author: Denis-Michael Lux	
# Date: 20. June 2015
#
# This software may be modified and distributed under the terms
# of the BSD license. See the LICENSE file for details.

# get gtest library
INCLUDE_DIRECTORIES( ${gtest_SOURCE_DIR}/include )

# collect unit test source files
SET(test_sources
	# files containing the actual tests
	${PROJECT_SOURCE_DIR}/test/SampleTest.cpp
	${PROJECT_SOURCE_DIR}/test/MatrixTest.cpp
)

# Build executable
ADD_EXECUTABLE        (  unit_tests ${test_sources}    )
TARGET_LINK_LIBRARIES (  unit_tests gtest_main UZLMath )