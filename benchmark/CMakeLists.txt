# Author: Denis-Michael Lux
# Date: 02. Juli 2015
#
# This software may be modified and distributed under the terms
# of the BSD license. See the LICENSE file for details.

# Creating binary and linking library
ADD_EXECUTABLE(         benchmark_soft_for ${PROJECT_SOURCE_DIR}/benchmark/benchmark_soft_for.cpp                       )
TARGET_LINK_LIBRARIES(  benchmark_soft_for UZLMath                                                                      )

ADD_EXECUTABLE(         benchmark_soft_for_speedup ${PROJECT_SOURCE_DIR}/benchmark/benchmark_soft_for_speedup.cpp       )
TARGET_LINK_LIBRARIES(  benchmark_soft_for_speedup UZLMath                                                              )

ADD_EXECUTABLE(         benchmark_soft_for_efficiency ${PROJECT_SOURCE_DIR}/benchmark/benchmark_soft_for_efficiency.cpp )
TARGET_LINK_LIBRARIES(  benchmark_soft_for_efficiency UZLMath                                                           )

ADD_EXECUTABLE(         benchmark_soft_inv ${PROJECT_SOURCE_DIR}/benchmark/benchmark_soft_inv.cpp                       )
TARGET_LINK_LIBRARIES(  benchmark_soft_inv UZLMath                                                                      )

ADD_EXECUTABLE(         benchmark_soft_inv_speedup ${PROJECT_SOURCE_DIR}/benchmark/benchmark_soft_inv_speedup.cpp       )
TARGET_LINK_LIBRARIES(  benchmark_soft_inv_speedup UZLMath                                                              )

ADD_EXECUTABLE(         benchmark_soft_inv_efficiency ${PROJECT_SOURCE_DIR}/benchmark/benchmark_soft_inv_efficiency.cpp )
TARGET_LINK_LIBRARIES(  benchmark_soft_inv_efficiency UZLMath                                                           )

ADD_EXECUTABLE(         benchmark_dwt_accuracy ${PROJECT_SOURCE_DIR}/benchmark/benchmark_dwt_accuracy.cpp               )
TARGET_LINK_LIBRARIES(  benchmark_dwt_accuracy UZLMath                                                                  )
