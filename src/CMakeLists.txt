cmake_minimum_required(VERSION 3.5)
project(CFD_LAB)

set(CMAKE_CXX_STANDARD 14)

FIND_PACKAGE(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

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
        LBMHelper.hpp
        LBMHelper.cpp
        ConfigParser.hpp
        VtkWriter.cpp
        VtkWriter.hpp
        freeSurface.cpp
        timeStep.hpp
        timeStep.cpp scenarios/Scenario.hpp scenarios/OnlyWater.cpp scenarios/OnlyWater.hpp scenarios/DamBreak.cpp scenarios/DamBreak.hpp scenarios/FallingDrop.cpp scenarios/FallingDrop.hpp scenarios/CornerDamBreak.cpp scenarios/CornerDamBreak.hpp)

add_executable(lbsim ${sources})
target_include_directories(lbsim PUBLIC ./)
target_link_libraries(lbsim m)
target_compile_options(lbsim PUBLIC -Wall -Wno-missing-braces -pedantic -Werror -g -O3 -march=native)
