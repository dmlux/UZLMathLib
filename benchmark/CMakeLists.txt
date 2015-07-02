# Author: Denis-Michael Lux
# Date: 11. January 2015
#
# This software may be modified and distributed under the terms
# of the BSD license. See the LICENSE file for details.

# Creating binary and linking library
ADD_EXECUTABLE( benchmark_soft_for ${PROJECT_SOURCE_DIR}/benchmark/benchmark_soft_for.cpp )
TARGET_LINK_LIBRARIES( benchmark_soft_for UZLMath )
