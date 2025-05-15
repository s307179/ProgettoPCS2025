cmake_minimum_required(VERSION 3.9)

project(PCSProject LANGUAGES CXX C VERSION 1.0.0)

# Set C++ Standard
################################################################################
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# IMPOSE WARNINGS ON DEBUG
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic-errors")

# Add dependencies
################################################################################

## Eigen3
find_package(Eigen3 CONFIG REQUIRED)
list(APPEND Polyhedron_LINKED_LIBRARIES PUBLIC Eigen3::Eigen)

## Threads
find_package(Threads REQUIRED)
list(APPEND Polyhedron_LINKED_LIBRARIES PRIVATE Threads::Threads)

## GTest
find_package(GTest REQUIRED)

list(APPEND Polyhedron_LINKED_LIBRARIES PRIVATE GTest::gtest)
list(APPEND Polyhedron_LINKED_LIBRARIES PRIVATE GTest::gtest_main)
list(APPEND Polyhedron_LINKED_LIBRARIES PRIVATE GTest::gmock)
list(APPEND Polyhedron_LINKED_LIBRARIES PRIVATE GTest::gmock_main)

# Insert Sources
################################################################################
add_subdirectory(src)

list(APPEND Polyhedron_SOURCES ${polyhedron_sources})
list(APPEND Polyhedron_HEADERS ${polyhedron_headers})
list(APPEND Polyhedron_INCLUDES ${polyhedron_includes})

add_subdirectory(src_test)
list(APPEND Polyhedron_SOURCES ${src_test_sources})
list(APPEND Polyhedron_HEADERS ${src_test_headers})
list(APPEND Polyhedron_INCLUDES ${src_test_includes})

add_subdirectory(ExportParaview)

list(APPEND Polyhedron_SOURCES ${src_paraview_sources})
list(APPEND Polyhedron_HEADERS ${src_paraview_headers})
list(APPEND Polyhedron_INCLUDES ${src_paraview_includes})

# Create executable
################################################################################
add_executable(${PROJECT_NAME}
      main.cpp
      ${Polyhedron_SOURCES}
      ${Polyhedron_HEADERS})


# To read our core data base
set(POLYHEDRON_TYPES "tetrahedron" "cube" "octahedron" "icosahedron" "dodecahedron")

foreach(poly_type ${POLYHEDRON_TYPES})
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Platonic_solids/${poly_type}.txt  
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endforeach()

target_link_libraries(${PROJECT_NAME} ${Polyhedron_LINKED_LIBRARIES})
target_include_directories(${PROJECT_NAME} PRIVATE ${Polyhedron_INCLUDES})
target_compile_options(${PROJECT_NAME} PUBLIC -fPIC)


add_executable(${PROJECT_NAME}_TEST
      main_test.cpp
      ${Polyhedron_SOURCES}
      ${Polyhedron_HEADERS}
      ${src_test_sources}
      ${src_test_headers})

target_link_libraries(${PROJECT_NAME}_TEST ${Polyhedron_LINKED_LIBRARIES})
target_include_directories(${PROJECT_NAME}_TEST PRIVATE ${Polyhedron_INCLUDES} ${src_test_includes})
target_compile_options(${PROJECT_NAME}_TEST PUBLIC -fPIC)