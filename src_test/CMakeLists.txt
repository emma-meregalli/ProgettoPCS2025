cmake_minimum_required(VERSION 3.10)
project(PCS_Project_Tests)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Percorsi alle librerie esterne (Eigen, GTest)
include_directories(/usr/include/eigen3)  

# Inclusione della directory dei test e dei sorgenti
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# File sorgenti del progetto (fuori dalla cartella test)
file(GLOB PROJECT_SOURCES
    ../*.cpp
)

# File sorgenti dei test (incluso main_test.cpp)
file(GLOB TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/main_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_utils/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_triangulation/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_dual/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_export/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_files/*.cpp
)

# Trova GoogleTest installato sul sistema
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Crea l'eseguibile dei test
add_executable(AllTests ${PROJECT_SOURCES} ${TEST_SOURCES})

# Collega solo gtest (non gtest_main!) + pthread
target_link_libraries(AllTests ${GTEST_LIBRARIES} pthread)

# Abilita i test
enable_testing()
add_test(NAME AllTests COMMAND AllTests)
