message("reading /test/CMakeLists.txt")

project(MatrixTests LANGUAGES C CXX)

# Specify the C++ standard to use:
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Specify the C standard to use:
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)

# Set the minimum version of cmake:
cmake_minimum_required(VERSION 3.0.0)

set(SOURCES main.cpp)

add_executable(matrix_tests ${SOURCES})

target_include_directories(matrix_tests PRIVATE ../src/include /usr/include/gtest)

target_link_libraries(matrix_tests PRIVATE matrix)
target_link_libraries(matrix_tests PRIVATE gtest)
target_link_libraries(matrix_tests PRIVATE pthread)
