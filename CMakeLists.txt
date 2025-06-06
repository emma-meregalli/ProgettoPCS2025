cmake_minimum_required(VERSION 3.9)

project(PolyhedralMesh LANGUAGES CXX C VERSION 1.0.0)

# Set C++ Standard
################################################################################
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# IMPOSE WARNINGS ON DEBUG
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic-errors")

# IMPOSE CXX FLAGS FOR WINDOWS
if (WIN32)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
endif (WIN32)

# Add dependencies
################################################################################

## Eigen3
find_package(Eigen3 CONFIG REQUIRED)
list(APPEND PolyhedralMesh_LINKED_LIBRARIES PUBLIC Eigen3::Eigen)

## Threads
find_package(Threads REQUIRED)
list(APPEND PolyhedralMesh_LINKED_LIBRARIES PRIVATE Threads::Threads)

## GTest
find_package(GTest REQUIRED)
list(APPEND PolyhedralMesh_LINKED_LIBRARIES PRIVATE GTest::gtest)
list(APPEND PolyhedralMesh_LINKED_LIBRARIES PRIVATE GTest::gtest_main)
list(APPEND PolyhedralMesh_LINKED_LIBRARIES PRIVATE GTest::gmock)
list(APPEND PolyhedralMesh_LINKED_LIBRARIES PRIVATE GTest::gmock_main)

# Insert Sources
################################################################################
file(GLOB_RECURSE polygonalMesh_sources CONFIGURE_DEPENDS
     src/Utils.cpp
     src/PolyhedralTriangulation.cpp)

file(GLOB_RECURSE polyhedralMesh_headers CONFIGURE_DEPENDS
     src/Utils.hpp
     src/PolyhedralTriangulation.hpp
     src/PolyhedralMesh.hpp)

include_directories(src)
list(APPEND PolyhedralMesh_SOURCES ${polygonalMesh_sources})
list(APPEND PolyhedralMesh_HEADERS ${polygonalMesh_headers})
list(APPEND PolyhedralMesh_INCLUDE src)

# Add ExportParaview sources
add_subdirectory(ExportParaview)

list(APPEND PolyhedralMesh_SOURCES ${src_paraview_sources})
list(APPEND PolyhedralMesh_HEADERS ${src_paraview_headers})
list(APPEND PolyhedralMesh_INCLUDE ${src_paraview_includes})

# Create main executable
################################################################################
add_executable(${PROJECT_NAME}
      main.cpp
      ${PolyhedralMesh_SOURCES}
      ${PolyhedralMesh_HEADERS})



# Target setup
################################################################################
target_link_libraries(${PROJECT_NAME} ${PolyhedralMesh_LINKED_LIBRARIES})
target_include_directories(${PROJECT_NAME} PRIVATE ${PolyhedralMesh_INCLUDE})
target_compile_options(${PROJECT_NAME} PUBLIC -fPIC)

# Create test executable
################################################################################
add_executable(${PROJECT_NAME}_TEST
      main_test.cpp
      ${PolyhedralMesh_SOURCES}
      ${PolyhedralMesh_HEADERS})

target_link_libraries(${PROJECT_NAME}_TEST ${PolyhedralMesh_LINKED_LIBRARIES})
target_include_directories(${PROJECT_NAME}_TEST PRIVATE ${PolyhedralMesh_INCLUDE})
target_compile_options(${PROJECT_NAME}_TEST PUBLIC -fPIC)
