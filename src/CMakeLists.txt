set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -g -O3 -march=native -fopenmp")
add_compile_options(-std=c++11)

add_library(weightedruns partitions.cxx pvalue.cxx)
target_link_libraries(weightedruns ${GSL_LIBRARIES})
target_include_directories(weightedruns PUBLIC "${GSL_INCLUDE_DIRS}")

add_executable(Partitions_Test partitions_TEST.cxx)
target_link_libraries(Partitions_Test weightedruns)
add_test(PartitionsTest Partitions_Test)

add_executable(pvalue_Test pvalue_TEST.cxx)
target_link_libraries(pvalue_Test weightedruns)
add_test(pvalueTest pvalue_Test)
# Local Variables:
# compile-command:"pushd ./build && cmake .. && make VERBOSE=1 && ctest --output-on-failure"
# End:
