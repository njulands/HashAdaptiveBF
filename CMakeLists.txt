cmake_minimum_required(VERSION 3.10)

# set the project name
project(CostSketch VERSION 0.1)

# add sub dir as lib
add_subdirectory(util)

# add the executable
add_executable(experiment experiment.cpp)

# link executable to library
target_link_libraries(experiment PUBLIC util)

# add the binary tree to the search path for include files
target_include_directories(experiment PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          "${PROJECT_SOURCE_DIR}/util"
                          )

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "-O3")