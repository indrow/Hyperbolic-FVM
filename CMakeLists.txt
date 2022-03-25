cmake_minimum_required(VERSION 3.16)
project(hertz)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -O3")

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

find_library(GMSH_LIB gmsh)

if (NOT GMSH_LIB)
    message(FATAL_ERROR "Could not find libgmsh")
endif ()

find_path(GMSH_INC gmsh.h)

if (NOT GMSH_INC)
    message(FATAL_ERROR "Could not find gmsh.h")
endif ()

include_directories(
        ${GMSH_INC}
        ${PROJECT_SOURCE_DIR}/src
)

add_executable(hertz src/main.cpp src/geo/CCoordinate.h src/geo/CVector.h src/geo/CLine.h src/common/CPi.h src/common/CType.h src/msh/CGmsh.h src/fvm/CBasis.h src/fvm/CCell.h src/fvm/CFace.h src/fvm/CGrid.h src/fvm/CNode.h src/fvm/CGridCreate.h src/solver/Scalar.h src/fvm/DCell.h src/fvm/DFace.h src/physics/Scalar.h src/schemes/RK.h)

target_link_libraries(hertz Eigen3::Eigen gmsh)