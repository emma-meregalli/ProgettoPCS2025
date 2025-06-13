cmake_minimum_required(VERSION 3.14)
project(PCS_Project_Tests)

set(CMAKE_CXX_STANDARD 17)

# === GoogleTest ===
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/release-1.12.1.zip
)
# Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

# === Eigen ===
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# === Include Directories ===
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
)

# === Source Files ===
set(SOURCES
    src/Utils.cpp
    src/PolyhedralTriangulation.cpp
    src/PolyhedralMesh.cpp
)

# === Test File ===
add_executable(runTests test/test.cpp ${SOURCES})

# === Link Libraries ===
target_link_libraries(runTests gtest_main Eigen3::Eigen)

# === Enable Testing ===
include(GoogleTest)
gtest_discover_tests(runTests)