cmake_minimum_required(VERSION 3.0)
project(runs)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(GSL REQUIRED)

enable_testing()

add_subdirectory(src)