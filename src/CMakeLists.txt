cmake_minimum_required(VERSION 3.5)
project(CFD_LAB)

set(CMAKE_CXX_STANDARD 14)

set(sources
        boundary.cpp
        boundary.hpp
        collision.cpp
        collision.hpp
        computeCellValues.cpp
        computeCellValues.hpp
        ConfigParser.hpp
        freeSurface.hpp
        initLB.cpp
        initLB.hpp
        LBDefinitions.hpp
        main.cpp
        streaming.cpp
        streaming.hpp
        parallel.cpp
        parallel.hpp
        LBMHelper.hpp
        LBMHelper.cpp
        ConfigParser.hpp
        VtkWriter.cpp
        VtkWriter.hpp freeSurface.cpp)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})

add_executable(lbsim ${sources})
target_include_directories(lbsim PUBLIC ./)
target_link_libraries(lbsim m ${MPI_CXX_LIBRARIES})
target_compile_options(lbsim PUBLIC -Wall -pedantic -Werror -g -O3 -march=native ${MPI_CXX_COMPILE_FLAGS})
