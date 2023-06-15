message("reading /CMakeLists.txt")

# Set the minimum version of cmake:
cmake_minimum_required(VERSION 3.0.0)

# Define the project:
set(CMAKE_PROJECT_NAME demo )
project(${CMAKE_PROJECT_NAME} VERSION 1.0 LANGUAGES C CXX)
message("CMAKE_PROJECT_NAME:  " ${CMAKE_PROJECT_NAME})
message("PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR})

# Specify the C standard to use:
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)

# Specify the C++ standard to use:
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Notify CMake of subdirectories:
add_subdirectory(src)
add_subdirectory(test)
