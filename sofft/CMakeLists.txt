# Author: Denis-Michael Lux
# Date: 22. Juli 2015
#
# This software may be modified and distributed under the terms
# of the BSD license. See the LICENSE file for details.

# Creating binary and linking library
ADD_EXECUTABLE(
    sofft_test
    ${PROJECT_SOURCE_DIR}/sofft/sofft_test.cpp
)
TARGET_LINK_LIBRARIES( sofft_test UZLMath )